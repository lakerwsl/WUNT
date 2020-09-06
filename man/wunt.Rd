\name{wunt}
\alias{wunt}
\title{Data-driven Uniform transformer.}
\usage{
wunt(x,z,y,type='constant',marginal=FALSE,extra=FALSE,xextra=NULL,method='Kernel',h=NULL,kernel='gaussian',degree=2,interact=TRUE)
}
\arguments{
\item{x}{A data matrix/vector}
\item{z}{A vector of treamtment assignment, with 1 for treated and 0 for control.}
\item{type}{Type of the smooth kernel function used to approximate the density of the covariates in the control group. type should be one of 'constant', 'linear' and 'quadratic'.}
\item{marginal}{An indicator for whether the Uniform transformer is constructed based on the marginal distribution or joint distribution. By defaul, marginal is FALSE, i.e., the Uniform transformer is constructed based on the joint distribution of x.}
\item{extra}{An indicator for whether there is extra unlabeled data to construct the Uniform transformer.}
\item{xextra}{Extra unlabeled data to construct the Uniform transformer. xextra must be specified if extra=T}
\item{method}{Method to estimate the density of the transformed data in the treated group. method should be 'Kernel' or 'Projection'}
\item{h}{The bandwidth to be used for the Kernel method. }
\item{kernel}{Type of kernel. kernel can be chosen from 'gaussian','gaussian4','gaussian6','gaussian8','gaussian10','epanechnikov','uniform','triangular','triweight','tricube','biweight',cosine','silverman'}
\item{degree}{Degree of polynomials for the basis functions of the Projection method.}
\item{interaction}{An indicator for whether interactions would be included for the basis functions.}
}
\value{
\item{Estimated average treatment effect on the treated group.}
\item{weight}{Constructed weights for the each subject in the control group.}
\item{control}{Transformed data for the control group. }
\item{treated}{Tranformed data for the treated group.}
}
\description{
This function constructs weights for each subject in the control group with a data-driven Uniform transformer, and estimates the average treatment effect on the treated group with the constructed weights.
}
\examples{
x1=rnorm(1000)
x2=rnorm(1000, mean=0.5)
x=cbind(x1,x2)
ix=sample(1:1000,300)
z=numeric(1000)
z[ix]=1
y=x1*x2
res1=wunt(x,z,y,type='constant',marginal=FALSE,method='Kernel')
res1$est
res2=wunt(x,z,y,type='constant',marginal=FALSE,method='Projection')
res2$est
}
