\name{getSamples}
\alias{getSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get samples }
\description{
  Fetch a subset of samples from a given data structure
}
\usage{
getSamples(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the matrix of sample data to subset. If \code{x} is a vector it is interpreted as a matrix. \code{x} can also be of class \code{"ExpressionSet"}. }
  \item{y}{ which sample(s) to use as data. Can be a vector of characters matching sample names, integers indicating which samples to choose, or a mixture of the two. If not supplied, all samples will be used. }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{getSamples-methods}}):
	\item{\code{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\code{order}}{ vector of characters, specifying on which column(s) to order the sample data. If \option{NULL} (default), the data will be returned without ordering enforced. }
	\item{\code{\dots}}{ other arguments passed are not handled at this time. }
    }
}
\value{
  Returns a matrix of values corresponding to a subset of samples from the data supplied, where columns correspond to samples.  Function halts if no samples to return.   
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{getSamples-methods}} }
\examples{
data(sample.ExpressionSet)
se.ABC <- getSamples(sample.ExpressionSet, c("A", "B", "C"), element="se.exprs")
se.ABC[1:10,]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
