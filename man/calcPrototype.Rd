\name{calcPrototype}
\alias{calcPrototype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate prototype }
\description{
  Calculates prototype (trimmed mean) across all samples
}
\usage{
calcPrototype(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric matrix, where each column represents a different sample }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{calcPrototype-methods}}):
    \describe{
	\item{\option{element}}{ which element of \code{AssayData} to use for a given \code{ExpressionSet} input (default is \option{"exprs"})}
	\item{\option{samples}}{ which samples to use as data. Can be a vector of characters matching sample names, integers indicating which samples to choose, or a mixture of the two. If \option{NULL} (default), all samples will be used.}
	\item{\option{center}}{ logical; if \option{TRUE} (default) samples will be mean-centered before protototype calculation. If \option{FALSE}, this mean-centering step will be skipped }
	\item{\option{trim}}{ the fraction (0 to 0.5, default is 0.1) of observations to be trimmed from each end of each row and column in \code{x} before the mean is computed. }
	\item{\option{verbose}}{ logical; if \option{TRUE} (default) progress will be output to screen. If \option{FALSE}, no output will be displayed.}
	\item{\option{\dots}}{ other arguments to be passed to \code{mean}. See \code{\link{mean}}.}
	}
    }
}
\value{
 Returns a vector of numerical data, representing the prototype ([trimmed] mean) of all samples in x.
}
\author{ Mark Reimers (\email{mreimers@vcu.edu}), Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{calcPrototype-methods}}, \code{\link{mean}} }
\examples{
#demo(pipeline,package="HELP")

x <- matrix(data=rep(1:1000,10),nrow=1000,ncol=10)
x <- x*(sample(1:100/100,size=10000,replace=TRUE))
x <- t(t(x)-1000*(1:10))
x[c(1:10,991:1000),]
x.avg <- calcPrototype(x)
x.avg[c(1:10,991:1000)]

#rm(x,x.avg)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
