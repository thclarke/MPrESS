% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpress.R
\name{plot.mpress}
\alias{plot.mpress}
\title{plot.mpress
makes a ggplot with the p-value and power for all the samples investigated an mpress object}
\usage{
\method{plot}{mpress}(x)
}
\description{
plot.mpress
makes a ggplot with the p-value and power for all the samples investigated an mpress object
}
\examples{
library(mpress);
#Loading in the  microbiome files
data(ChinaData);
data(SpainData);

#Use summary() to examine the data loaded

#Use plot() to see the plot of the power and p-value in the different sample numbers
plot(spain.ibs.power)
}
