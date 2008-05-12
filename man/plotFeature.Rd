\name{plotFeature}
\alias{plotFeature}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot feature versus two-color intensity  }
\description{
  Graphical display of featureData (ex: fragment size) versus two-color signal intensity data
}
\usage{
plotFeature(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of numerical data to be plotted, with two columns (one for each signal channel). \code{x} can also be of class \code{"ExpressionSet"}. }
  \item{y}{ an additional vector of numerical data to be used for feature. If \code{y} is missing, the function will attempt to fill a value from \code{featureData} in \code{x}. }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{plotFeature-methods}}):
	\item{\code{element.x}}{ which element of \code{AssayData} to use (for signal channel 1) for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\code{element.y}}{ which element of \code{AssayData} to use (for signal channel 2) for a given \code{ExpressionSet} input (default is \option{"exprs2"}) }
	\item{\code{sample}}{ which element of \code{sampleNames} to use as data (default is 1). Can be a character matching a sample name or simply an integer indicating which sample to choose. }
	\item{\code{feature}}{ which element of \code{featureData} to use as plotting feature (default is 1). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. }
	\item{\code{feature.random}}{ which element of \code{featureData} to use to identify random probes (default is \option{"TYPE"}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. }
	\item{\code{which.random}}{ an integer vector specifying which rows of data correspond to random probes. if \option{NULL} (default), the function will attempt to identify random probes using \code{featureData}. }
	\item{\code{random.flag}}{ a character specifying the label for random probes in \code{feature.random} from \code{featureData}. Default is \option{"RAND"}. }
	\item{\code{na.rm}}{ logical; if \option{TRUE} (default), missing values are removed from \code{x}. If \option{FALSE} any missing values cause an error.}
	\item{\code{limit}}{ numerical input specifying the maximum number of points to plot (default is 10,000). if \option{NULL}, all points will be used.}
	\item{\code{cutoff}}{ a numerical input specifying the value below which signal intensities from channel 1 can be considered "failed" probes. If \option{NULL} (default), the function will attempt to calculate a cutoff from random probe information.}
	\item{\code{cutoff2}}{ a numerical input specifying the value below which signal intensities from channel 2 can be considered "failed" probes. If \option{NULL} (default), the function will attempt to calculate a cutoff from random probe information.}
	\item{\code{main}}{ an overall title for the plot: see \code{\link{title}}. }
	\item{\code{xlab}}{ a title for the x axis (default is \option{"Fragment size (bp)"}): see \code{\link{title}}. }
	\item{\code{ylab}}{ a title for the y axis for signal channel 1 (default is \option{"log(MspI)"}): see \code{\link{title}}. }
	\item{\code{ylab2}}{ a title for the y axis for signal channel 2 (default is \option{"log(HpaII)"}): see \code{\link{title}}. }
	\item{\code{cex}}{ numerical value (default is 0.2) giving the amount by which plotting text and symbols should be scaled relative to the default.  }
	\item{\code{\dots}}{ other arguments to be passed to \code{plot}. See \code{\link{plot}}.}  
 }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{plotFeature-methods}} }
\examples{
#demo(pipeline,package="HELP")

msp1 <- sample(8000:16000/1000, size=1000)
msp1 <- msp1[order(msp1)]
hpa2 <- sample(8000:16000/1000, size=1000)
hpa2 <- hpa2[order(hpa2)]
size <- sample((1:1000)*1.8+200, size=1000)
rand <- which.min(abs(msp1-quantile(msp1, 0.25)))
plotFeature(cbind(msp1, hpa2), size, which.random=(rand-20):(rand+20), main="Random")

#rm(msp1, hpa2, size, rand)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
