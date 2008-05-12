\name{calcPrototype-methods}
\docType{methods}
\alias{calcPrototype-methods}
\alias{calcPrototype,missing-method}
\alias{calcPrototype,ExpressionSet-method}
\alias{calcPrototype,vector-method}
\title{ Calculate prototype (methods)}
\description{
  Methods for calculating prototype ([trimmed] mean) across all samples
}
\section{Methods}{
\describe{

\item{x = "missing"}{ Handle empty function call }

\item{x = "ExpressionSet"}{ Handle input of an object of class \code{ExpressionSet}. Derive data from \code{AssayData}. }

\item{x = "vector"}{ Handle vector input as a matrix }
}}
\author{ Mark Reimers (\email{mreimers@vcu.edu}), Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{calcPrototype}} }
\keyword{methods}
