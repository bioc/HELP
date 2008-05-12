\name{combineData-methods}
\docType{methods}
\alias{combineData-methods}
\alias{combineData,missing,missing,missing-method}
\alias{combineData,vector,missing,missing-method}
\alias{combineData,vector,missing,vector-method}
\alias{combineData,vector,vector,missing-method}
\alias{combineData,vector,vector,vector-method}
\alias{combineData,ExpressionSet,missing,missing-method}
\alias{combineData,ExpressionSet,vector,missing-method}
\title{ Combine data (methods)}
\description{
  Methods for calculating trimmed and/or weighted means of groups of rows in a given data matrix.
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing", w = "missing"}{ Handle empty function call }

\item{x = "vector", y = "missing", w = "missing"}{ Handle partially empty function call. Reinterpret with default parameters instead of missing values. }

\item{x = "vector", y = "missing", w = "vector"}{ Handle partially empty function call. Reinterpret with default parameters instead of missing values.  }

\item{x = "vector", y = "vector", w = "missing"}{ Handle partially empty function call. Reinterpret with default parameters instead of missing values.  }

\item{x = "vector", y = "vector", w = "vector"}{ Handle input of three vectors specifying data, grouping, and weighting information, respectively. Note that the data and weighting inputs are handled as matrices.  }

\item{x = "ExpressionSet", y = "missing", w = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive grouping and weighting data from \code{featureData} and \code{AssayDataElement}, respectively. }

\item{x = "ExpressionSet", y = "vector", w = "missing"}{ Handle input of an object of class \code{ExpressionSet} and a vector specifying grouping information. Derive weighting data from code{AssayDataElement}. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{combineData}}}
\keyword{methods}
