\name{exprs2}
\alias{exprs2}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{exprs2<-}
\title{ Retrieve microarray data (for signal channel 2) from ExpressionSets }
\description{
  Access (and/or assign) data for signal channel 2 in a given ExpressionSet object
}
\usage{
exprs2(object)
exprs2(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ Object of class \code{ExpessionSet} }
  \item{value}{ Matrix with rows representing features and columns representing samples }
}
\value{
  \code{exprs2} returns a (usually large!) matrix of values
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\examples{
#demo(pipeline,package="HELP")
}
% \seealso{ \code{\link{exprs}}, \code{\link{ExpressionSet-class}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
