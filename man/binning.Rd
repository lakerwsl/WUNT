\name{binning}
\alias{binning}
\title{Dividing a vector of numbers uniformly into bins.}
\usage{
binning(v,nbin,range=c(0,1))
}
\arguments{
\item{v}{A vector of numbers.}
\item{nbin}{Number of bins.}
\item{range}{A vector of two numbers indicating the range of elements in v.}
}
\value{
\item{b}{A vector with length length(v), giving which bin each element of v belongs to.}
\item{s}{A vector with length nbin, giving the lower bound defining each bin.}
\item{e}{A vector with length nbin, giving the upper bound defining each bin.}
}
\description{
This function divides v uniformly into nbin and finds the lower bound and upper bound for each bin.

Of limited interest to most users; function binning() would typically be called by some other functions.
}
\examples{
binning(1:10,2,c(0,10))
}
