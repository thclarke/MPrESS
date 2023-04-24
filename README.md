# # MPrESS

MPrESS (Microbiome Power Estimation with Sampling and Simulation) is an R-package which can estimate the power necessary to calculate the significant difference between two sample types in a microbiome as represented by a phyloseq class object.

## Installing MPrESS

MPrESS can be installed using the install_git command in devtools.  The example data is also pulled for demonstration purposes.
```
R
>library(devtools)
>library(git2r)
>install_git("https://github.com/JCVenterInstitute/MPrESS.git")
>quit()
```

if this is not working, it is also the possible to download mpress_1.0.0.tar.gz from github and install MPrESS from the zip file.

## Using MPrESS

Using the R-console to use the MPrESS function to calculates the number of samples to detect a difference between microbiomes with two different metadata values. Below we will show the default version of running it with the example data from Chavdra et al 2016 and .

### Using MPrESS in the R console

MPrESS is centered around two primary functions to load and analyze the data, with three primary output variables that allow for analysis and reporting of the clustering. For more detailed descriptions of all the GGRaSP R functions, please examine the vigenettes.

To start with, simply load the library and Spain IBS and China Geographic 
>library(mpress);
>#Loading in the Chinese metadata file
>data(ChinaData);
>data(SpainData);
>
>china.power = power.est(china.full, "State", c("Yunnan", "Guangxi Zhuang"));
>#Use summary() to examine the data loaded
>summary(china.power)
>
>#Other microbiomes: china.trim from data(ChinaData): unnamed OTUs removed
>#Other microbiomes: spain.ibs.full from data(SpainData): all OTUs removed
>#Other microbiomes: spain.ibs.trim from data(SpainData): unnamed OTUs removed
>
>ibs.power = power.est(spain.ibs.trim, "Disease", c("H", "M"));
>
>Use plot() to see the plot of the power and p-value in the different sample numbers
>plot(ibs.power)

##Citing MPrESS

When using MPrESS in your analysis, please cite MPrESS using the github page.
