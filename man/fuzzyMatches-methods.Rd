\name{fuzzyMatches-methods}
\docType{methods}
\alias{fuzzyMatches-methods}
\alias{fuzzyMatches,missing,missing-method}
\alias{fuzzyMatches,vector,missing-method}
\alias{fuzzyMatches,vector,NULL-method}
\alias{fuzzyMatches,vector,vector-method}
\alias{fuzzyMatches,NULL,vector-method}
\title{ Fuzzy matching (methods)}
\description{
 Methods for matching and reinterpreting a vector in terms of a second vector, essentially using the second vector as a key to interpreting the first.
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "vector", y = "missing"}{ Handle empty function call }

\item{x = "vector", y = "NULL"}{ Handle empty function call }

\item{x = "vector", y = "vector"}{  Handle input of two vectors. }

\item{x = "NULL", y = "vector"}{ Handle empty function call }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{fuzzyMatches}} }
\keyword{methods}
