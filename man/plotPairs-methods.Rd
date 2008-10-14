\name{plotPairs-methods}
\docType{methods}
\alias{plotPairs-methods}
\alias{plotPairs,missing-method}
\alias{plotPairs,matrix-method}
\alias{plotPairs,ExpressionSet-method}
\title{ Plot tree-pairs (methods)}
\description{
 Methods for pairwise comparison of samples producing a matrix of scatterplots and a corresponding dendrogram
}
\section{Methods}{
\describe{

\item{x = "missing"}{ Handle empty function call }

\item{x = "matrix"}{ Handle matrix input }

\item{x = "ExpressionSet"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from \code{AssayData}. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{plotPairs}} }
\keyword{methods}
