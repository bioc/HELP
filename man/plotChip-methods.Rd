\name{plotChip-methods}
\docType{methods}
\alias{plotChip-methods}
\alias{plotChip,missing,missing,missing-method}
\alias{plotChip,matrix,missing,missing-method}
\alias{plotChip,ExpressionSet,missing,missing-method}
\alias{plotChip,ExpressionSet,vector,missing-method}
\alias{plotChip,ExpressionSet,ExpressionSet,missing-method}
\alias{plotChip,vector,vector,ExpressionSet-method}
\alias{plotChip,vector,vector,vector-method}
\title{ Plot chip image (methods)}
\description{
 Methods for graphic display of spatially-linked data, particularly applicable for microarrays
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing", z = "missing"}{ Handle empty function call }

\item{x = "matrix", y = "missing", z = "missing"}{ Handle matrix input, extract information, and reinterpret function call with appropriate vectors }

\item{x = "ExpressionSet", y = "missing", z = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive both data and position information from a single object. }

\item{x = "ExpressionSet", y = "vector", z = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive position information from this object, but the corresponding data from vector input. }

\item{x = "ExpressionSet", y = "ExpressionSet", z = "missing"}{ Handle input of two objects of class \code{ExpressionSet}. Derive position information and data from each one, respectively. }

\item{x = "vector", y = "vector", z = "ExpressionSet"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from this object, but the corresponding position information from vector input. }

\item{x = "vector", y = "vector", z = "vector"}{ Handle input of three vectors. Derive X and Y positions and data from each one, respectively. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}), Mark Reimers (\email{mreimers@vcu.edu}) }
\seealso{ \code{\link{plotChip}} }
\keyword{methods}
