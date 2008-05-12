\name{quantileNormalize-methods}
\docType{methods}
\alias{quantileNormalize-methods}
\alias{quantileNormalize,missing,missing-method}
\alias{quantileNormalize,matrix,missing-method}
\alias{quantileNormalize,vector,missing-method}
\alias{quantileNormalize,vector,ExpressionSet-method}
\alias{quantileNormalize,vector,vector-method}
\alias{quantileNormalize,ExpressionSet,missing-method}
\alias{quantileNormalize,ExpressionSet,vector-method}
\alias{quantileNormalize,ExpressionSet,ExpressionSet-method}
\title{ Quantile normalization (methods)}
\description{
 Methods for applying quantile normalization to multiple bins of data, divided by a sliding window approach
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "matrix", y = "missing"}{ Handle matrix input, reinterpret function call with two vector input if matrix has two columns, otherwise handle as empty function call }

\item{x = "vector", y = "missing"}{ Handle empty function call }

\item{x = "vector", y = "ExpressionSet"}{ Handle input of an object of class \code{ExpressionSet}. Derive binning information from this class but use data from a vector input. }

\item{x = "vector", y = "vector"}{ Handle input of two vectors specifying data and binning information, respectively. }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive both data and binning information from a single object. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from this class but use binning information from a vector input. }

\item{x = "ExpressionSet", y = "ExpressionSet"}{ Handle input of two objects of class \code{ExpressionSet}. Derive data and binning information from each one, respectively. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{quantileNormalize}} }
\keyword{methods}
