\name{createWiggle-methods}
\docType{methods}
\alias{createWiggle-methods}
\alias{createWiggle,missing,missing-method}
\alias{createWiggle,ExpressionSet,missing-method}
\alias{createWiggle,ExpressionSet,vector-method}
\alias{createWiggle,vector,vector-method}
\title{ Create wiggle track (methods)}
\description{
 Methods for creating wiggle tracks
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive features from \code{FeatureData}. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Derive features from vector input }

\item{x = "vector", y = "vector"}{ Handle vector input }
}}
\author{Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{createWiggle}} }
\keyword{methods}
