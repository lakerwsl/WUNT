\name{transformer}
\alias{transformer}
\title{Data-driven Uniform transformer.}
\usage{
transformer(x,newx,marginal=FALSE,type='constant')
}
\arguments{
\item{x}{A data matrix/vector to construct the Uniform transformer}
\item{newx}{A new data matrix/vector to be transformed.}
\item{marginal}{An indicator for whether the Uniform transformer is constructed based on the marginal distribution or joint distribution. By defaul, marginal is FALSE, i.e., the Uniform transformer is constructed based on the joint distribution of x.}
\item{type}{Type of the smooth kernel function used to approximate the density of the covariates in the control group. type should be one of 'constant', 'linear' and 'quadratic'.}
}
\value{
\item{xt}{Transformed data of x with the constructed Uniform transformer.}
\item{newxt}{Transformed data of newx with the constructed Uniform transformer.}
}
\description{
This function constructs a data-driven Uniform transformer to transform the distribution of x to a Uniform distribution, and returns the transformed data for the original data x and a new data newx.
}
\examples{
res=transformer(rnorm(100),rnorm(100,mean=0.5))
res$xt
res$newxt
}
