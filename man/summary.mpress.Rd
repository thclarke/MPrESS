% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpress.R
\name{summary.mpress}
\alias{summary.mpress}
\title{summary.mpress
writes the power value and the estimation type to the screen}
\usage{
\method{summary}{mpress}(x)
}
\description{
summary.mpress
writes the power value and the estimation type to the screen
}
\examples{
library(mpress);
#Loading in the  microbiome files
data(ChinaData);
data(SpainData);

summary(china.trim.power)
}
