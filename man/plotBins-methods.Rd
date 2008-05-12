\name{plotBins-methods}
\docType{methods}
\alias{plotBins-methods}
\alias{plotBins,missing,missing-method}
\alias{plotBins,matrix,missing-method}
\alias{plotBins,vector,missing-method}
\alias{plotBins,vector,ExpressionSet-method}
\alias{plotBins,vector,vector-method}
\alias{plotBins,ExpressionSet,missing-method}
\alias{plotBins,ExpressionSet,vector-method}
\alias{plotBins,ExpressionSet,ExpressionSet-method}
\title{ Plot bins (methods) }
\description{
 Methods for plotting densities of multiple bins of data, divided by a sliding window approach
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "matrix", y = "missing"}{ Handle matrix input, reinterpret function call with two vector input if matrix has two columns, otherwise handle as empty function call }

\item{x = "vector", y = "missing"}{ Handle empty function call }

\item{x = "vector", y = "ExpressionSet"}{ Handle input of an object of class \code{ExpressionSet}. Derive binning information from this class but use data from a vector input. }

\item{x = "vector", y = "vector"}{ Handle input of two vectors specifying data and binning information, respectively. }

\item{x = "matrix", y = "matrix"}{ Handle matrix input, reinterpret function call with vector input }

\item{x = "matrix", y = "vector"}{ Handle matrix input, reinterpret function call with vector input }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive both data and binning information from a single object. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from this class but use binning information from a vector input. }

\item{x = "ExpressionSet", y = "ExpressionSet"}{ Handle input of two objects of class \code{ExpressionSet}. Derive data and binning information from each one, respectively. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{plotBins}} }
\keyword{methods}
