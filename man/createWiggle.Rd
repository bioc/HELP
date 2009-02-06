\name{createWiggle}
\alias{createWiggle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create wiggle track }
\description{
  Create and write a wiggle track (UCSC Genome Browser format) to flat file
}
\usage{
createWiggle(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ matrix of numerical data, where each column represents data for an individual wiggle track. \code{x} can also be of class \code{"ExpressionSet"}.  }
  \item{y}{  an additional matrix of numerical data with columns corresponding to chr, start, and end, respectively. }
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{createWiggle-methods}}):
     \describe{
	\item{\option{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\option{feature.chr}}{ which element of \code{featureData} to use as chromosomal information (default is \option{"CHR"}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. }
	\item{\option{feature.start}}{ which element of \code{featureData} to use as start positions (default is \option{"START"}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. }
	\item{\option{feature.stop}}{ which element of \code{featureData} to use as end positions (default is \option{"STOP"}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. }
	\item{\option{samples}}{ which sample(s) to use as data. Can be a vector of characters matching sample names, integers indicating which samples to choose, or a mixture of the two. If \option{NULL} (default), all samples will be used. }
	\item{\option{colors}}{ vector of colors, indicates which colors to use for which wiggle track }
	\item{\option{file}}{ location of file to write wiggle track information; if \option{""}, wiggle track prints to the standard output connection: see \code{\link{cat}}. }
	\item{\option{append}}{ logical; if \option{TRUE}, the output is appended to an existent wiggle track file. If \option{FALSE} (default), a new file with a new header is created and any existing file of the same name is destroyed. }
	\item{\option{na.rm}}{ logical; if \option{TRUE} (default), missing values are removed from data. If \option{FALSE} any missing values cause an error }
	\item{\option{sep}}{ a string used to separate columns. Using \option{sep = "\t"} (default) gives tab-delimited output. }
	\item{\option{\dots}}{ other arguments to be passed to \code{cat}. See \code{\link{cat}}.}
    }
  }
}
\references{ UCSC Genome Browser, \url{http://genome.ucsc.edu/goldenPath/help/customTrack.html}:
Kent, W.J., Sugnet, C. W., Furey, T. S., Roskin, K.M., Pringle, T. H., Zahler, A. M., and Haussler, D. The Human Genome Browser at UCSC. \emph{Genome Res}. \bold{12}(6), 996-1006 (2002). }
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{write}}, \code{\link{cat}} }
\examples{
#demo(pipeline,package="HELP")

chr <- rep("chr1", 500)
start <- (1:500)*200
end <- start+199
data <- sample(5*(1:10000/10000)-2, size=500)
data <- cbind(data, abs(data), -1*data)
colnames(data) <- c("data", "abs", "opposite")
createWiggle(data, cbind(chr, start, end))

#rm(chr, start, end, data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }
