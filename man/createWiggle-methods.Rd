\name{createWiggle-methods}
\docType{methods}
\alias{createWiggle-methods}
\alias{createWiggle,missing,missing-method}
\alias{createWiggle,ExpressionSet,missing-method}
\alias{createWiggle,ExpressionSet,matrix-method}
\alias{createWiggle,vector,matrix-method}
\alias{createWiggle,matrix,matrix-method}
\title{ Create wiggle track (methods)}
\description{
 Methods for creating wiggle tracks
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive features from \code{FeatureData}. }

\item{x = "ExpressionSet", y = "matrix"}{ Handle input of an object of class \code{ExpressionSet}. Derive features from matrix input }

\item{x = "vector", y = "matrix"}{ Handle vector input }

\item{x = "matrix", y = "matrix"}{ Handle matrix input }
}}
\author{Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{createWiggle}} }
\keyword{methods}
