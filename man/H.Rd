\name{H}
\alias{H}
\title{Integtal of the smooth kernel function used to approximate the density of the covariates in the control group.}
\usage{
H(x,type='constant')
}
\arguments{
\item{x}{A real number.}
\item{type}{Type of the smooth kernel function used to approximate the density of the covariates in the control group. type should be one of 'constant', 'linear' and 'quadratic'.}
}
\value{
Integtal of the smooth kernel function used to approximate the density of the covariates in the control group, evaluted at x.
}
\description{
This function calculates the intergal of the smooth kernel function used to approximate the density of the covariates in the control group, at x.

Of limited interest to most users; function H() would typically be called by some other functions.
}
\examples{
H(0.1)
}
