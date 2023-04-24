#' @title power.calc
#' @description power.calc() calculates the number of samples to dectect a difference between microbiomes with two different metadata values.
#'
#' @param in.phyloseq Required. The collection of microbiome samples to use. Required to be a class phyloseq data
#' @param metadata.var Required. The column name in the metadata values that is being compared for the power.
#' @param metadata.vals Required. Vector with the values to be compared.
#' @param start. The sample size at which to start the tests. Default is 5.
#' @param alpha Alpha value used to test the microbiome differences. Default is 0.05 
#' @param beta Beta value used to test proportion of estimated samples that contain significant differences. Default is 0.95
#' @param dist.metric String giving the distance metric to use between the microbiome samples. Default "bray" (Bray-Curtis). Available options can be seen with distanceMethodList. 
#' @param deseq.val Value to use for the DESeq2 run to trim the taxa. if the value is 0,use all taxa (Default); if the value is between 0 and 1, all taxa with an FDR below the value will be selected; otherwise, if the value is >=1, that number of top taxa by FDR will be selected
#' @param n.rep Maximum number of replicates to run before final check if above the beta threshold. Default = 100.
#' @param burn.in Value to burning in for the binomial test. Default is 10
#' @param verbose Flag to prints intermediate values to the screen. Default = F
#' @param binom Flag to use the binomial test to stop low mapping values. Default = F
#' @param switch.val sample number at which to automatically swap between sampling and simulation. Default = 0, or only swap once sample number exceeds the number of samples
#' @param seed.val Value to set seed value for troubleshooting. Default = 0 (not set)
#' @return Returns a class MPrESS variable with the sample number in addition to runtime information for use in verification
#' @examples #The following data is from Chavda et al 2016 which phylotyped Enterobacter genomes
#'
#' #Additional printing and plotting options are availible with plot() and print(). For more information refer to ?plot.mpress and ?print.mpress
#'@examples #The following data is from Zhang et al 2015 using Microbiome data from multiple sites in China
#' # Our example uses the data underpinning the tree shown in Figure 2
#' 
#' library(mpress);
#' #Loading in the Chinese metadata file
#' data(ChinaData);
#' data(SpainData);
#'
#` china.power = power.est(china.full, "State", c("Yunnan", "Guangxi Zhuang"));
#'#Use summary() to examine the data loaded
#' summary(china.power)
#'
#' #Other microbiomes: china.trim from data(ChinaData): unnamed OTUs removed
#' #                  spain.ibs.trim from
#'
#' #Use plot() to see the plot of the power and p-value in the different sample numbers
#' plot(china.power)
#' @export

power.est = function(in.phyloseq, metadata.var, metadata.vals, start.n=5, alpha = 0.05, beta = 0.95, dist.metric="bray", deseq.val = 0, n.rep=100, burn.in = 10, verbose = T, binom = F, switch.val = 0, seed.val=0, tax.val ="genus")
{
  if (seed.val != 0)
  {
    set.seed(seed.val)
  }
  #Set test type to permanova as currently the only option. Maybe be changed to a input variable if more tests are added later
  test.type="permanova";
  #Step one, verify variables and that metadata.var and metadata.vals are in the 
  if (is.null(phyloseq) | class(in.phyloseq) != "phyloseq")
  {
    cat("Must include phyloseq-class variable as in.phyloseq")
    return(NULL);
  }
  if (is.null(metadata.var) | !(metadata.var %in% colnames(sample_data(in.phyloseq))))
  {
    cat(paste("Cannot find metadata variable ", metadata.var," in in.phyloseq. Quitting..."));
        return(NULL);
  }
  if (length(metadata.vals) > 2)
  {
    cat("Only two metadata variable values can be compared...")
    return(NULL)
  }
  
  if (length(metadata.vals) < 2)
  {
    cat("Two metadata variable values are  required...")
    return(NULL)
  }

  #Step one and a ha;f
  if (!tax.val %in% rank_names(in.phylosq)){
    cat(paste("Cannot locate ", tax.val, " in in.phyloseq. Using the full OTU table in submitted data\n", sep=""))
  }
  else{
      in.phyloseq = tax_glom(in.phyloseq, tax.val)
  }
  #Step two, get minimal counts for each metadata.val
  min.cnt = .get_min_metadata_counts(in.phyloseq, metadata.var, metadata.vals);
  if (min.cnt <= 0)
  {
    cat(paste("Cannot find metadata values ", metadata.var," in in.phyloseq. Quitting..."))
        return(NULL);
  }
  #Step three: trim by deseq2 if needed
  if (deseq.val > 0)
  {
    in.phyloseq = .trim_taxa_by_deseq(in.phyloseq, metadata.var, metadata.vals, deseq.val, verbose);
  }
  #Step three get the beta value for the first sample number
  cur.num = start.n;
  run.num = 0;
  p.val = c();
  shape = NA;
  alpha.list = c();
  est.type = "sample";
  if (verbose==TRUE)
  {
    cat(paste("# of Samples", "# of Replicates", "Estimation Type", "Mean p-value",  "Power", "# of Signifcant Replicates", "\n", sep=","))
  }
  ret.df = data.frame(stringsAsFactors=F);
  while(length(p.val) == 0 | mean(p.val < alpha) < beta)
  {
    alpha.list = c();
    p.val = c()
    run.num = 0;
    p.v = 1;
    
    while(run.num <= n.rep & p.v > 0.001 & sum(p.val > alpha) < (1-beta) * n.rep)
    {
      if ((cur.num > min.cnt | (cur.num > switch.val & switch.val != 0) ) & est.type =="sample")
      {
        est.type = "simulation";
        shape = .get_sim_shapes(in.phyloseq, metadata.var, metadata.vals);
      }
      if (est.type == "sample")
      {
        tmp.phyloseq = .get_sampled_otus(in.phyloseq, k = cur.num, metadata.var, metadata.vals);
      }
	  if (est.type == "simulation")
      {
        tmp.phyloseq = .get_simulated_otus(in.phyloseq, metadata.var, metadata.vals, cur.num, shape)
      }
      alpha.list = c(alpha.list, estimate_richness(tmp.phyloseq)$Shannon);
      dist.val = .get_distance_value(tmp.phyloseq, tolower(dist.metric))
      test.1 <- .get_test_pvalue(tmp.phyloseq, metadata.var, dist.val, test.type);
      run.num = run.num + 1;
      p.val = c(p.val, test.1);
      if (run.num >= burn.in & binom ==TRUE)
      {
        p.v = binom.test(sum(p.val < alpha), run.num, beta)$p.value;
      }
      
    }
    if (verbose==TRUE)
    {
      cat(paste(cur.num, run.num, est.type, mean(p.val),  mean(p.val < alpha), sum(p.val < alpha), "\n"))
    }
    ret.df = rbind(ret.df, data.frame(num.samples = cur.num, num.rep = run.num, estimation.type = est.type, p.val.mean = mean(p.val), p.val.sd = sd(p.val), power.val = mean(p.val < alpha), alpha.mean = mean(alpha.list), alpha.sd = sd(alpha.list)))
    cur.num = cur.num + 1;
  }
  
  power.ret = new("mpress", sample.number= cur.num-1, test = test.type, distance=dist.metric, estimation.type=est.type, in.data = in.phyloseq, runtime.info = ret.df);
  
  return(power.ret);
} 

#' An S4 class representing the MPrESS data and output
#'
#' @slot sample.number Integer given the number of runs to get the power
#' @slot test String given the test type
#' @slot distance String given the distance metrix
#' @slot estimation.type String given the estimation type
#' @slot in.data A phyloseq object used to create the simulation
#' @slot runtime.info A data.frame with the stats for each run

setClass("mpress", representation(sample.number="numeric", test = "character", distance="character", estimation.type="character", in.data = "phyloseq", runtime.info ="data.frame"));

setMethod("show", "mpress", function(object)summary(object))

#' print.mpress
#' makes a table with the runtime information for the samples investigated an mpress object
#' 
#' @examples
#' library(mpress);
#' #Loading in the microbiome files
#' data(ChinaData);
#' data(SpainData);
#'
#` spain.ibs.power = power.est(spain.ibs.full, "Disease`, c("H", "M"));
#'#Use summary() to examine the data loaded
#'
#' #Use print(...) to see the table of the run info in the different sample numbers
#' print(spain.ibs.power)
#' #Use print(..., y="otu") to see the table of the run info in the different sample numbers
#' print(spain.ibs.power, y="otu")

#' @export
print.mpress <- function(x, type="runtime")
{
  if (y == "runtime"){
    write.table(x@runtime.info, sep="\t", quote=F);
  }
  else {
    if (y == "otu"){
      cnt.table <- data.frame(otu_table(x@in.data), stringsAsFactors=F);
      write.table(cnt.table, sep="\t", quote=F);
    }
    else{
      write("Type is not recognized. Please use either runtime or otu\n");
    }
}

#' summary.mpress
#' writes the power value and the estimation type to the screen
#' 
#' @examples
#' library(mpress);
#' #Loading in the  microbiome files
#' data(ChinaData);
#' data(SpainData);
#'
#` china.trim.power = power.est(china.trim, "State",c("Yunnan", "Guangxi Zhuang"));
# '#Use summary() to examine the data loaded
#' summary(china.trim.power)
#' @export
summary.mpress <- function(x)
{
  cat(paste("", x@sample.number, " using ", x@estimation.type, "\n", sep=""));
}

#' plot.mpress
#' makes a ggplot with the p-value and power for all the samples investigated an mpress object
#' 
#' @examples
#' library(mpress);
#' #Loading in the  microbiome files
#' data(ChinaData);
#' data(SpainData);
#'
#` spain.ibs.power = power.est(spain.ibs.full, "Disease`, c("H", "M"));
#'#Use summary() to examine the data loaded
#'
#' #Use plot() to see the plot of the power and p-value in the different sample numbers
#' plot(spain.ibs.power)
#' @export
plot.mpress <- function(x)
{
  g.1 <- ggplot(data=x@runtime.info, aes(x = num.samples, y = p.val.mean, group=estimation.type, linetype = estimation.type, color="Mean P-Value"))+geom_smooth() + geom_smooth(aes(x = num.samples, y = power.val, color="Power", group=estimation.type, linetype = estimation.type))+xlab("Sample Number")+ylab("P-value/Power")
  g.1 <- g.1 + geom_point(data=x@runtime.info, aes(x = num.samples, y = power.val, shape=estimation.type))
  return(g.1)
}


.trim_taxa_by_deseq = function(test.phylo, metadata.var, metadata.vals, alpha, verbose.in)
{
  x1 = c();
  if (verbose.in == TRUE)
  {
    cat(paste("Starting DESeq...\n", sep=""))
  }
  test.phylo.2 = prune_samples(sample_names(test.phylo)[(sapply(sample_data(test.phylo)[,metadata.var], function(x){ x %in% metadata.vals}))], test.phylo);
  cnt.table <- data.frame(otu_table(test.phylo.2), stringsAsFactors=F);
  tax.table <- data.frame(tax_table(test.phylo.2), stringsAsFactors=F);
  metadata.table <- data.frame(var = sample_data(test.phylo.2)[,metadata.var]);
  new.metadata.vals = gsub(pattern = "([^0-9A-Za-z\\_\\.]+)", "", metadata.vals)
  metadata.table[,metadata.var] = gsub(pattern = "([^0-9A-Za-z\\_\\.]+)", "", metadata.table[,metadata.var])
  
  tmp.dds <- DESeqDataSetFromMatrix(cnt.table, colData=metadata.table, design = as.formula(paste("~", metadata.var)));
  tmp.dds.run <- DESeq(tmp.dds, fitType = "local", quiet= T)
  tmp.dds.res <- results(tmp.dds.run);
  hit = 0;
  max.pval = -1;
  if (alpha > 0 & alpha < 1)
  {
    if (sum(tmp.dds.res$padj < alpha & !is.na(tmp.dds.res$padj)) >1)
    {
      hit = 1;
      tmp.list = rownames(tmp.dds.res)[tmp.dds.res$padj < alpha & !is.na(tmp.dds.res$padj)];
      max.pval = max(tmp.dds.res$padj[tmp.dds.res$padj < alpha & !is.na(tmp.dds.res$padj)]);
      test.phylo.2 <- prune_taxa(taxa_names(test.phylo) %in% tmp.list, test.phylo);
    }
  }
  if (alpha >= 1)
  {
    hit = 1;
    tmp.list = rownames(tmp.dds.res)[order(tmp.dds.res$padj, na.last=T)][1:alpha];
    max.pval = sort(tmp.dds.res$padj, na.last=T)[alpha];
    test.phylo.2 <- prune_taxa(taxa_names(test.phylo) %in% tmp.list, test.phylo);
  }
  if (verbose.in == TRUE)
  {
    cat(paste("DESeq Finished...\nTrimmed down to ", length(tmp.list), " taxa with a max FDR of ", max.pval,"\n\n\n", sep=""))
  }
  return(test.phylo.2);
}


.get_sim_shapes = function(in.phylo, metadata.var, metadata.vals)
{
  shp.list = list();
  shp.list[["name"]] = list();
  shp.list[["shape"]] = list();
  for (st in metadata.vals)
  {
    tmp.in <- t(as.data.frame(otu_table(in.phylo))[, sample_data(in.phylo)[,metadata.var]==st])
    shp.list[["name"]][[st]] = colnames(tmp.in)[colSums(tmp.in)>0]
    shp.list[["shape"]][[st]] <- dirmult(tmp.in)$gamma
  }
  return(shp.list);
}

.get_simulated_otus = function(in.phylo, metadata.var, metadata.vals, k, shape.list)
{
  test.phylo = NULL
  samp.data <- data.frame(Var = c(), stringsAsFactors = F);
  samp.num = list()
  for (st in metadata.vals)
  {
    shp <- shape.list[["shape"]][[st]];
    samp.num[[st]] = sample(sample_sums(in.phylo), k, replace = T);
    otu.samp <- Dirichlet.multinomial(samp.num[[st]] , shp)  
    colnames(otu.samp) = shape.list[["name"]][[st]]
    rownames(otu.samp) = paste(st, 1:k, sep="");
    samp.data.tmp <- data.frame( Var=rep(st,k), stringsAsFactors = F);
    rownames(samp.data.tmp) = paste(st, 1:k, sep="");
    
    samp.data = rbind(samp.data, samp.data.tmp);
    if (!class(test.phylo) == "phyloseq")
    {
      test.phylo = phyloseq(otu_table(t(otu.samp), taxa_are_rows = T), tax_table(tax_table(in.phylo)[taxa_names(in.phylo) %in% colnames(otu.samp),]))  
    }
    else
    {
      tmp.phylo = phyloseq(otu_table(t(otu.samp), taxa_are_rows = T), tax_table(tax_table(in.phylo)[taxa_names(in.phylo) %in% colnames(otu.samp),]))            
      test.phylo = merge_phyloseq(test.phylo, tmp.phylo)
    }
  }
  colnames(samp.data) = c(metadata.var);
  sample_data(test.phylo) <- samp.data;
  if (!is.null( access(in.phylo, "phy_tree")))
  {
    phy_tree(test.phylo) = phy_tree(in.phylo);
  }
  return(test.phylo);
}

.get_sampled_otus = function(in.phylo, metadata.var, metadata.vals, k)
{
  test.phylo = NULL
  n.1 = c();
  for (st in metadata.vals)
  {
    n.1 <- c(n.1, sample(sample_names(in.phylo)[sample_data(in.phylo)[,metadata.var] == st],k, replace = F))
  }
  test.phylo = prune_samples(sample_names(in.phylo) %in% n.1, in.phylo)
  return(test.phylo)
}



.get_distance_value = function(test.phylo, dist.type)
{
  d.list = phyloseq::distanceMethodList;
  dist.val = NA;
  test.phylo = phyloseq::transform_sample_counts(test.phylo, function(OTU) OTU /sum(OTU) )
  if (dist.type %in% d.list$vegdist)
  {
    dist.val = vegan::vegdist(t(phyloseq::otu_table(test.phylo)), dist.type);
  }
  else
  {	
    dist.val = phyloseq::distance(test.phylo, dist.type);
  }
  return(dist.val);
}

.get_test_pvalue = function(test.phylo, metadata.var, dist.val, test.type)
{
  p.val = NA;
  if (test.type=="permanova")
  {
    test.val = adonis(as.formula(paste("dist.val ~ ", metadata.var, sep="")),  as(sample_data(test.phylo), "data.frame"));
    p.val = test.val$aov.tab$`Pr(>F)`[1];
  }
  return(p.val);
}

.get_min_metadata_counts = function(in.phylo, metadata.var, metadata.vals)
{
  min.cnt = min(sapply(metadata.vals, function(st) { sum(sample_data(in.phylo)[,metadata.var] == st)}));
  return(min.cnt);
}