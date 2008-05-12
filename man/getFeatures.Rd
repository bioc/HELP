\name{getFeatures}
\alias{getFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get features (methods) }
\description{
  Fetch a subset of features from a given data structure
}
\usage{
getFeatures(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the matrix of feature data to subset. If \code{x} is a vector it is interpreted as a matrix. \code{x} can also be of class \code{"ExpressionSet"} or an \code{"AnnotatedDataFrame"} object. }
  \item{y}{ which feature(s) to use. Can be a vector of characters matching feature names, integers indicating which features to choose, or a mixture of the two. If not supplied (or if equivalent to \option{"*"}), all features will be used. }
  \item{\dots}{ other arguments passed are not handled at this time. }
}
\value{
  Returns a matrix of values corresponding to a subset of features from the data structure supplied, where columns correspond to features.  Function halts if no features to return.   
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{getFeatures-methods}} }
\examples{
data(sample.ExpressionSet)
df <- data.frame(x=1:500,y=501:1000, row.names=featureNames(sample.ExpressionSet))
featureData(sample.ExpressionSet) <- new("AnnotatedDataFrame", data=df, dimLabels=c("featureNames", ""))
getFeatures(sample.ExpressionSet, "y")[1:10]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }
