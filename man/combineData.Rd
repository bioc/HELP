\name{combineData}
\alias{combineData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Combine data }
\description{
  Calculate trimmed and/or weighted means of groups of rows in a given data matrix.
}
\usage{
combineData(x, y, w, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric matrix containing the values whose trimmed and/or weighted mean is to be computed. Each column is treated independently. }
  \item{y}{ a vector describing the discrete groups used to divide the elements of x.  If y is missing then all elements of x are handled together. }
  \item{w}{ a matrix of weights the same dimensions as x giving the weights to use for each element of x. If w is missing then all elements of x are given the same weight. }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{combineData-methods}}): 
	\item{\code{trim}}{the fraction (0 to 0.5, default is 0) of observations to be trimmed from each group of rows in \code{x} according to \code{y}. }
	\item{\code{na.rm}}{logical; if \option{TRUE}, missing values are removed from \code{x} and \code{y} and \code{z}. If \option{FALSE} any missing values cause an error.}
	\item{\code{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\code{feature.group}}{ which element of \code{featureData} to use as binning variable (default is \code{NULL}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. See \code{\link{getFeatures}}. }
	\item{\code{element.weight}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input.  If \code{NULL} (default), weighting is not performed. }
	\item{\code{feature.weight}}{ which element of \code{featureData} to use as weighting variable (default is \code{NULL}). Can be a character matching \code{varLabel} or simply an integer indicating which feature to choose. See \code{\link{getFeatures}}. }
	\item{\code{samples}}{ which samples to use as data. Can be a vector of characters matching sample names, integers indicating which samples to choose, or a mixture of the two. If \option{NULL} (default), all samples will be used.}
	\item{\code{\dots}}{ other arguments not handled at this time.}  
  }
}
\value{
  Returns a matrix of combined numerical data, where each row represents the summary of a group of elements from the corresponding column in x.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\note{ Each column in data matrix treated separately.}
\seealso{ \code{\link{combineData-methods}}, \code{\link{mean}}, \code{\link{weighted.mean}} }
\examples{
#demo(pipeline,package="HELP")

x <- 1:100
combineData(x,w=x/100)
weighted.mean(x,w=x/100)

y <- sample(c("a","b","c",1:3),size=100,replace=TRUE)
combineData(cbind(x,x,2*x),y,trim=0.5)

#rm(x,y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith }
