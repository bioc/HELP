\name{plotChip}
\alias{plotChip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot chip image }
\description{
  Graphic display of spatially-linked data, particularly applicable for microarrays
}
\usage{
plotChip(x, y, z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of numerical data determining x-coordinates of data on chip. \code{x} can also handle \code{ExpressionSet} (see \code{\link{plotChip-methods}} for more parameter details). }
  \item{y}{ vector of numerical data determining y-coordinates of data on chip }
  \item{z}{ the vector of numerical data to be plotted }
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{plotChip-methods}}):
     \describe{
	\item{\code{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\code{sample}}{ which element of \code{sampleNames} to use as
	data (default is 1). Can be a character matching a sample name or
	simply an integer indicating which sample to choose.} 
	\item{\code{feature.x}}{ which element of \code{featureData} to use as
	X coordinate (default is \code{"X"}). Can be a character matching
	\code{varLabel} or simply an integer indicating which feature to
	choose. } 
	\item{\code{feature.y}}{ which element of \code{featureData} to use as
	Y coordinate (default is \code{"Y"}). Can be a character matching
	\code{varLabel} or simply an integer indicating which feature to
	choose. } 
	\item{\code{na.rm}}{ logical; if \option{TRUE}, missing values are
	removed from \code{x}, \code{y}, and \code{z}. If \option{FALSE}
	(default) any missing values cause an error. } 
	\item{\code{main}}{ an overall title for the plot: see \code{\link{title}}. }
	\item{\code{xlab}}{ a title for the x axis: see \code{\link{title}}. }
	\item{\code{ylab}}{ a title for the y axis: see \code{\link{title}}. }
	\item{\code{colors}}{ vector of colors specifying the color scheme to
	use (default is \code{rev(rainbow(n=20, start=0, end=1/3))}). Also
	determines the resolution of \code{z} such that the more colors that
	are used allow finer discrimination of differences in \code{z}. } 
	\item{\code{range}}{ vector of numerical data of length 2 (default is
	\code{c(NA, NA)}) specifying range used to color-code data in \code{z}
	} 
	\item{\code{nrows}}{ numerical input specifying the number of rows by
	which to divide the chip; default is \option{NULL} which skips the
	division of data into blocks and results in individual spot resolution
	} 
	\item{\code{ncols}}{ numerical input specifying the number of columns
	by which to divide the chip; default is \option{NULL} which skips the
	division of data into blocks and results in individual spot resolution
	} 
	\item{\code{\dots}}{ other arguments to be passed to \code{plot}. See \code{\link{plot}}. }
	}
    }
}
\author{ Reid F. Thompson \email{rthompso@aecom.yu.edu}, Mark Reimers \email{mreimers@vcu.edu} }
\seealso{ \code{\link{plotChip-methods}} }
\examples{
#demo(pipeline,package="HELP")

x <- rep(1:100,100)
y <- rep(1:100,each=100)
z <- x*(1001:11000/1000)
z <- z-mean(z)
z <- z*(sample(1:10000/10000)+1)
plotChip(x,y,z,main="Curved gradient",xlab="x",ylab="y")

plotChip(x,y,sample(1:10000,size=10000),colors=gray(0:50/50),range=c(1,10000),main="Random noise")

#rm(x,y,z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
