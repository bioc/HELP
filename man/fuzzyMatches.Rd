\name{fuzzyMatches}
\alias{fuzzyMatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fuzzy matching }
\description{
  Match and reinterpret a vector in terms of a second vector, essentially using the second vector as a key to interpret the first.
}
\usage{
fuzzyMatches(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector, the values to be matched. }
  \item{y}{ vector, the values to be matched against. }
  \item{\dots}{ Arguments to be passed to methods (see
  \code{\link{getSamples-methods}}):
    \describe{
	\item{\code{strict}}{ logical indicating whether or not to allow restrict matching. If \option{FALSE}, numerical indices in \code{x} can reference values in \code{y}. If \option{TRUE} (default), only exact matches of values in \code{x} and \code{y} will be used. }
	\item{\code{keep}}{ logical indicating whether or not to preserve non-matching values from \code{x} (when \code{strict = FALSE}). If \option{TRUE} (default), all values in \code{x} will be returned with those that match values in \code{y} replaced by the corresponding values in \code{y}. If \option{FALSE}, non-matching values will be removed.  }
	\item{\code{alias}}{ logical indicating whether or not to return aliased values (default is \option{TRUE}). If \option{FALSE}, numerical indices will be returned with a value of \code{nomatch} for non-matching values in \code{x}. }
	\item{\code{match.all}}{ character value, specifying a special value to be interpreted as a match for ALL values in \code{y} (default is \option{"*"}). Any occurence of \code{match.all} in \code{x} will be replaced by all values in \code{y}. }
	\item{\code{nomatch}}{ numerical, specifying a value for non-matching values in code{x} when \code{strict = FALSE}, \code{keep = TRUE}, and \code{alias = FALSE}.}
	\item{\code{na.rm}}{ a logical value indicating whether NA values should be stripped before the computation proceeds (default is \option{TRUE}). If \option{FALSE} any missing values in \code{x} will cause an error and missing values in code{y} may cause unexpected behavior. }
	\item{\code{\dots}}{ other arguments not handled at this time.  }
    }
  }
}
\details{
  This function employs multiple stages of matching between two vectors.  First, the values in \code{x} are matched against \code{y} to find any exact matches. Next, numerical values in \code{x} are used to retrieve the corresponding positions in \code{y}. All unmatched values may be retained or dropped (depending on the value of \code{keep}), and a list of unique values is returned. Note that a value of \code{match.all} in \code{x} will be interpreted as a match for ALL values in \code{y}, and therefore replaced with the contents of \code{y}.
}
\value{
	Returns a vector of unique values in \code{x}, that match values in \code{y} according to the parameters described above.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{fuzzyMatches-methods}}, \code{\link{match}} }
\examples{
a <- c(1, "four", "missing")
b <- c("one", "two", "three", "four")
fuzzyMatches(a, b)
fuzzyMatches(a, b, strict=FALSE)
fuzzyMatches(a, b, strict=FALSE, alias=FALSE)
fuzzyMatches(a, b, strict=FALSE, keep=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
