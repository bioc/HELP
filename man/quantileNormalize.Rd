\name{quantileNormalize}
\alias{quantileNormalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile normalization }
\description{
  Apply quantile normalization to multiple bins of data, divided by a sliding window approach
}
\usage{
quantileNormalize(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the vector of numerical data to be normalized. If \code{x} is a matrix it is interpreted as a vector. \code{x} can also be of class \code{"ExpressionSet"}. }
  \item{y}{ an additional vector of numerical data to be used for binning. If \code{y} is a matrix it is interpreted as a vector. \code{y} can also be of class \code{"ExpressionSet"}. }
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{quantileNormalize-methods}}):
     \describe{
	\item{\code{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\code{sample}}{ which element of \code{sampleNames} to use as data (default is 1). Can be a character matching a sample name or simply an integer indicating which sample to choose. See \code{\link{getSamples}}. }
	\item{\code{feature}}{ which element of \code{featureData} to use as binning variable (default is 1). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. See \code{\link{getFeatures}}. }
	\item{\code{num.bins}}{ number of bins (default is 10) used to divide the data }
	\item{\code{num.steps}}{ number of steps (default is 3) used to create bin offsets, resulting in bins of sliding windows }
	\item{\code{mode}}{ the binning mode to be used. This must be either \option{"continuous"} (default) or \option{"discrete"}. \option{"continuous"} mode will divide the data into density-dependent bins.  \option{"discrete"} mode will divide the data uniformly by binning data values. }
	\item{\code{type}}{ an integer between 1 and 9 (default is 7) selecting one of the nine quantile algorithms: see \code{\link{quantile}}. }
	\item{\code{na.rm}}{ logical; if \option{TRUE}, missing values are removed from x and y. If \option{FALSE} any missing values cause an error. }
	\item{\code{\dots}}{ other arguments to be passed to \code{quantile}. See \code{\link{quantile}}.}
	}
   }
}
\value{
  Returns a vector of normalized numerical data according to input parameters.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{quantileNormalize-methods}}, \code{\link{quantile}} }
\examples{
#demo(pipeline,package="HELP")

x <- rep(1:100,10)+10*rep(1:10,each=100)
y <- rep(1:20,each=50)
d <- density(quantileNormalize(x,y,num.bins=20,num.steps=1,mode="discrete"))
plot(density(x))
lines(d$x,d$y/3,col="red")

#rm(x,y,d)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
