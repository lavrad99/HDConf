% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDConf.R
\name{Ridge_pvalue}
\alias{Ridge_pvalue}
\title{Function for extracting p-values from Ridge regression using the method described in Cule, Vineis and De Iorio (2011)}
\usage{
Ridge_pvalue(x, y, conf, method = "Cule-DeIorio", n_lambdas = 200, t_test = T)
}
\arguments{
\item{x}{is the exposure vector}

\item{y}{is the response vector}

\item{conf}{is the confound matrix}

\item{method}{refers to method for selecting ridge penalty ("Cule-DeIorio" (default) or "CV")}

\item{n_lambdas}{refers is the number of lambdas used in cross-validation (only used when method = "CV")}

\item{t_test}{is a boolean determining whether to use t-tests (default) or z-tests}
}
\value{
the p-value and the sign of the coefficient for x
}
\description{
Function for extracting p-values from Ridge regression using the method described in Cule, Vineis and De Iorio (2011)
}
