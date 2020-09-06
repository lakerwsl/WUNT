\name{Kernel}
\alias{Kernel}
\title{Kernel density estimation.}
\usage{
Kernel(x,h,kernel='gaussian')
}
\arguments{
\item{x}{Data matrix or data vector to calculate the kernel estimate.}
\item{h}{The bandwidth to be used. }
\item{kernel}{Type of kernel. kernel can be chosen from 'gaussian','gaussian4','gaussian6','gaussian8','gaussian10','epanechnikov','uniform','triangular','triweight','tricube','biweight',cosine','silverman'}
}
\value{
A vector of the estimated density values.
}
\description{
This function calculates the kernel density estimates of x with kernel choice kernel and bandwidth h. If the number of covariates is large, higher order kernels are recommended.

Of limited interest to most users; function Kernel() would typically be called by some other functions.
}
\examples{
Kernel(x=cbind(runif(10),rnorm(10)),h=1/10,kernel='gaussian')
}
