\name{readPairs-methods}
\docType{methods}
\alias{readPairs-methods}
\alias{readPairs,missing,missing,missing-method}
\alias{readPairs,vector,missing,missing-method}
\alias{readPairs,vector,vector,missing-method}
\alias{readPairs,character,character,ExpressionSet-method}
\alias{readPairs,character,character,character-method}
\title{ Read Nimblegen .pair files (methods)}
\description{
 Methods for extracting data from corresponding files in the Nimblegen .pair format.
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing", z = "missing"}{ Handle empty function call }

\item{x = "vector", y = "missing", z = "missing"}{ Handle single vector input. If two values specified in vector, reinterpret function call with two character inputs. Otherwise, handle as empty function call. }

\item{x = "vector", y = "vector", z = "missing"}{ Handle two vector input.  If vectors of unit length, reinterpret function call with two character inputs. Otherwise, handle as improper function call. }

\item{x = "character", y = "character", z = "ExpressionSet"}{ Handle two character vector inputs, each specifiying a filename to use when reading pair information. Pair data will be written to an \code{ExpressionSet} object. }

\item{x = "character", y = "character", z = "character"}{ Handle two character vector inputs, each specifiying a filename to use when reading pair information. Pair data will be written to a database. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{readPairs}} }
\keyword{methods}
