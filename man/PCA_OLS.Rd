% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDConf.R
\name{PCA_OLS}
\alias{PCA_OLS}
\title{Function for high-dimensional "OLS-PCA"}
\usage{
PCA_OLS(x, y, conf, prop_VE = 0.95, min_DF = 20)
}
\arguments{
\item{x}{is the exposure vector}

\item{y}{is the response vector}

\item{conf}{is the confound matrix}

\item{prop_VE}{is the \% of Variance Explained criterion for including PCs}

\item{min_DF}{is the minimal number of degrees of freedom in the OLS-PCA approach}
}
\value{
the p-value and the sign of the coefficient for x
}
\description{
Function for high-dimensional "OLS-PCA"
}
\details{
If the VE criterion yields degrees of freedom below min_DF, the function defaults to HD_OLS
}
