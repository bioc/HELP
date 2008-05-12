\name{readDesign-methods}
\docType{methods}
\alias{readDesign-methods}
\alias{readDesign,missing,missing,missing-method}
\alias{readDesign,vector,missing,missing-method}
\alias{readDesign,vector,vector,missing-method}
\alias{readDesign,character,character,ExpressionSet-method}
\alias{readDesign,character,character,character-method}
\title{ Read NimbleGen design files (methods)}
\description{
 Methods for extracting array design information from corresponding files in the Nimblegen .ndf and .ngd formats.
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing", z = "missing"}{ Handle empty function call }

\item{x = "vector", y = "missing", z = "missing"}{ Handle single vector input. If two values specified in vector, reinterpret function call with two character inputs. Otherwise, handle as empty function call. }

\item{x = "vector", y = "vector", z = "missing"}{Handle two vector input.  If vectors of unit length, reinterpret function call with two character inputs. Otherwise, handle as improper function call. }

\item{x = "character", y = "character", z = "ExpressionSet"}{ Handle two character vector inputs, each specifiying a filename to use when reading design information. Design information will be written to an \code{ExpressionSet}.}

\item{x = "character", y = "character", z = "character"}{Handle two character vector inputs, each specifiying a filename to use when reading design information. Design information will be written to a database.}
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{readDesign}} }
\keyword{methods}
