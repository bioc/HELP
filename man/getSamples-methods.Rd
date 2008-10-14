\name{getSamples-methods}
\docType{methods}
\alias{getSamples-methods}
\alias{getSamples,missing,missing-method}
\alias{getSamples,ExpressionSet,missing-method}
\alias{getSamples,ExpressionSet,NULL-method}
\alias{getSamples,ExpressionSet,vector-method}
\alias{getSamples,vector,missing-method}
\alias{getSamples,vector,NULL-method}
\alias{getSamples,vector,vector-method}
\alias{getSamples,matrix,missing-method}
\alias{getSamples,matrix,NULL-method}
\alias{getSamples,matrix,vector-method}
\title{ Get samples (methods)}
\description{
 Methods for fetching subsets of samples from various data structures
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Select data for all samples. }

\item{x = "ExpressionSet", y = "NULL"}{ Handle input of an object of class \code{ExpressionSet}. Select data for all samples. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Select data for a subset of samples. }

\item{x = "vector", y = "missing"}{ Handle input of a vector (interpreted as a matrix). Select data for all samples. }

\item{x = "vector", y = "NULL"}{ Handle input of a vector (interpreted as a matrix). Select data for all samples. }

\item{x = "vector", y = "vector"}{ Handle input of two vectors specifying data and sample subset information, respectively. }

\item{x = "matrix", y = "missing"}{ Handle input of a matrix. Select data for all samples. }

\item{x = "matrix", y = "NULL"}{ Handle input of a matrix. Select data for all samples. }

\item{x = "matrix", y = "vector"}{ Handle input of a matrix and a vector specifying data and sample subset information, respectively. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{getSamples}} }
\keyword{methods}
