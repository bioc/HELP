\name{plotBins}
\alias{plotBins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot bins }
\description{
  Plot densities of multiple bins of data, divided by a sliding window approach
}
\usage{
plotBins(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the vector of numerical data to be plotted. If \code{x} is a matrix it is interpreted as a vector. \code{x} can also be of class \code{"ExpressionSet"}. }
  \item{y}{ an additional vector of numerical data to be used for binning. If \code{y} is a matrix it is interpreted as a vector. \code{y} can also be of class \code{"ExpressionSet"}.}
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{plotBins-methods}}):
     \describe{
	\item{\option{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\option{sample}}{ which element of \code{sampleNames} to use as data (default is 1). Can be a character matching a sample name or simply an integer indicating which sample to choose. See \code{\link{getSamples}}. }
	\item{\option{feature}}{ which element of \code{featureData} to use as binning variable (default is 1). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. See \code{\link{getFeatures}}.}
	\item{\option{num.bins}}{ number of bins (default is 10) used to divide the data }
	\item{\option{num.steps}}{ number of steps (default is 3) used to create bin offsets, resulting in bins of sliding windows }
	\item{\option{mode}}{ the binning mode to be used. This must be either \option{"continuous"} or \option{"discrete"}. \option{"continuous"} mode will divide the data into density-dependent bins.  \option{"discrete"} mode will divide the data uniformly by binning data values. }
	\item{\option{show.avg}}{ logical; if \option{TRUE}, plots overall density in addition to densities per bin. If \option{FALSE} (default), overall density plot is omitted. }
	\item{\option{main}}{ an overall title for the plot: see \code{\link{title}}. }
	\item{\option{xlab}}{ a title for the x axis: see \code{\link{title}}. }
	\item{\option{ylab}}{ a title for the y axis: see \code{\link{title}}. }
	\item{\option{na.rm}}{ logical; if \option{TRUE} (default), missing values are removed from x and y. If \option{FALSE} any missing values cause an error.}
	\item{\option{\dots}}{ other arguments to be passed to \code{plot}. See \code{\link{plot}}.}
	}
   }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{plotBins-methods}}, \code{\link{density}}, \code{\link{quantile}} }
\examples{
#demo(pipeline,package="HELP")

x <- 1:1000
y <- sample(1:50,size=1000,replace=TRUE)
plotBins(x,y,show.avg=TRUE,main="Random binning data",xlab="1:1000")

#rm(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
