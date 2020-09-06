\name{binning_search}
\alias{binning_search}
\title{Finding which bin each element of a vector belongs to.}
\usage{
binning_search(v,start,end)
}
\arguments{
\item{v}{A vector of numbers.}
\item{start}{lower bound of each bin}
\item{end}{upper bound of each bin}
}
\value{
A vector with length length(v), giving which bin each element of v belongs to.
}
\description{
This function determines which bin each element of belongs to based on the lower bound and upper bound of each bin. The range of the ith bin is (start[i],end[i]].

Of limited interest to most users; function binning_search() would typically be called by some other functions.
}
\examples{
binning_search(1:10,c(0,5),c(5,10))
}
