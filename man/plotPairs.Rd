\name{plotPairs}
\alias{plotPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot tree-pairs }
\description{
  Pairwise comparison of samples producing a matrix of scatterplots and a corresponding dendrogram
}
\usage{
plotPairs(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric matrix, where each column represents a different sample. \code{x} can also be of class \code{"ExpressionSet"}. }
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{plotPairs-methods}}):
     \describe{
	\item{\option{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"}) }
	\item{\option{samples}}{ which samples to use as data. Can be a vector of characters matching sample names, integers indicating which samples to choose, or a mixture of the two. If \option{NULL} (default), all samples will be used.}
	\item{\option{scale}}{ logical value indicating whether sample branch lengths should be scaled by distance (default is \option{TRUE}) }
	\item{\option{groups}}{ logical value indicating whether the samples should be organized and color-coded by group (default is \option{TRUE}) }
	\item{\option{dist.method}}{ the distance measure to be used. This must be one of \option{"euclidean"} (default), \option{"maximum"}, \option{"manhattan"}, \option{"canberra"}, \option{"binary"} or \option{"minkowski"}. Any unambiguous substring can be given: see \code{\link{dist}} for more details. }
	\item{\option{hclust.method}}{ the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of \option{"ward"} (default), \option{"single"}, \option{"complete"}, \option{"average"}, \option{"mcquitty"}, \option{"median"} or \option{"centroid"}: see \code{\link{hclust}} for more details. }
	\item{\option{k}}{ an integer scalar or vector with the desired number of groups. If \option{NULL} (default), grouping will rely instead on distance measurements: see \code{\link{cutree}} for more details. }
	\item{\option{\dots}}{ other arguments to be passed to \code{pairs} or \code{dist}. See \code{\link{pairs}}, \code{\link{dist}}.}
  	}
   }
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{plotPairs-methods}}, \code{\link{dist}}, \code{\link{hclust}}, \code{\link{dendrogram}}, \code{\link{cutree}}, \code{\link{pairs}} }
\examples{
#demo(pipeline,package="HELP")

x <- sample(1:10000,size=10000)
x <- cbind(x,x+5,x*sample((1000:2000)/1000,size=10000,replace=TRUE),sample(-1*(1:10000),size=10000))
colnames(x) <- c("x","x+5","spread","random")
plotPairs(x)

#rm(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
