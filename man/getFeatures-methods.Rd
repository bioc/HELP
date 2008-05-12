\name{getFeatures-methods}
\docType{methods}
\alias{getFeatures-methods}
\alias{getFeatures,missing,missing-method}
\alias{getFeatures,ExpressionSet,missing-method}
\alias{getFeatures,ExpressionSet,NULL-method}
\alias{getFeatures,ExpressionSet,vector-method}
\alias{getFeatures,AnnotatedDataFrame,missing-method}
\alias{getFeatures,AnnotatedDataFrame,NULL-method}
\alias{getFeatures,AnnotatedDataFrame,vector-method}
\alias{getFeatures,vector,missing-method}
\alias{getFeatures,vector,NULL-method}
\alias{getFeatures,vector,vector-method}
\title{ Get features (methods)}
\description{
 Methods for fetching a subset of features from a given data structure
}
\section{Methods}{
\describe{

\item{x = "missing", y = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet", y = "missing"}{ Handle input of an object of class \code{ExpressionSet}. Select all feature data. }

\item{x = "ExpressionSet", y = "NULL"}{ Handle input of an object of class \code{ExpressionSet}. Select all feature data. }

\item{x = "ExpressionSet", y = "vector"}{ Handle input of an object of class \code{ExpressionSet}. Select a subset of features. }

\item{x = "AnnotatedDataFrame", y = "missing"}{ Handle input of an \code{AnnotatedDataFrame} object. Select all feature data. }

\item{x = "AnnotatedDataFrame", y = "NULL"}{ Handle input of an \code{AnnotatedDataFrame} object. Select all feature data. }

\item{x = "AnnotatedDataFrame", y = "vector"}{ Handle input of an \code{AnnotatedDataFrame} object. Select a subset of features. }

\item{x = "vector", y = "missing"}{ Handle input of a vector (interpreted as a matrix). Select all feature data }

\item{x = "vector", y = "NULL"}{ Handle input of a vector (interpreted as a matrix). Select all feature data }

\item{x = "vector", y = "vector"}{ Handle input of two vectors specifying feature data and feature subset information, respectively. }
}}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{getFeatures}} }
\keyword{methods}
