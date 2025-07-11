% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDConf.R
\name{Naive_ElasticNet}
\alias{Naive_ElasticNet}
\title{Function performing high-dimensional estimation via Naive ElasticNet selection}
\usage{
Naive_ElasticNet(
  x,
  y,
  conf,
  criterion = "min",
  n_lambdas = 200,
  n_alphas = 10,
  min_DF = 20
)
}
\arguments{
\item{x}{is the exposure vector}

\item{y}{is the response vector}

\item{conf}{is the confound matrix}

\item{criterion}{refers to the criterion for choosing lambda, "min" (default) or "1se"}

\item{n_lambdas}{is the number of lambdas tested for selection using 5-fold cross-validation}

\item{n_alphas}{is the number of alphas tested for selection using 5-fold cross-validation}

\item{min_DF}{is the minimal number of degrees of freedom in the final OLS model}
}
\value{
the p-value and the sign of the coefficient for x
}
\description{
Function performing high-dimensional estimation via Naive ElasticNet selection
}
\details{
If the full model has degrees of freedom above min_DF, the function uses HD_OLS at the last step
}
