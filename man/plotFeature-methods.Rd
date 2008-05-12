\name{plotFeature-methods}
\docType{methods}
\alias{plotFeature-methods}
\alias{plotFeature,missing,missing-method}
\alias{plotFeature,ExpressionSet,missing-method}
\alias{plotFeature,ExpressionSet,vector-method}
\alias{plotFeature,matrix,vector-method}
\title{ Plot feature versus two-color intensity (methods)}
\description{
 Methods for plotting featureData (ex: fragment size) versus two-color signal intensity data
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Derive both data and feature information from a single object. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from this class but use feature values from a vector input. }

\item{x = "matrix", y = "vector"}{ Handle matrix input, where each of two columns in matrix represents data from one signal channel. Feature data is derived from values from a vector input. }
}}
\keyword{methods}
