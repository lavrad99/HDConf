% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDConf.R
\name{pinv_SVD}
\alias{pinv_SVD}
\title{Function for performing the Moore-Penrose pseudo-inverse using SVD}
\usage{
pinv_SVD(X, numerical_zero = 1e-07)
}
\arguments{
\item{X}{is a matrix}

\item{numerical_zero}{is the level at which an eigenvalue is considered zero}
}
\value{
The Moore-Penrose pseudo-inverse of X
}
\description{
Function for performing the Moore-Penrose pseudo-inverse using SVD
}
