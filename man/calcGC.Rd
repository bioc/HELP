\name{calcGC}
\alias{calcGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate GC percent }
\description{
  Function to calculate GC percent from a nucleotide sequence input
}
\usage{
calcGC(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ characters containing nucleotide sequence (ex: \option{"ATCGGAA"}) or an object of class \code{ExpressionSet} }
  \item{\dots}{ Other arguments passed to methods:
	\code{\item{allow}}{ vector of characters specifying what other characters to allow in sequence (default is \option{"N"}) }
	}
}
\value{
  Returns a numerical value (from \option{0} to \option{1}) indicating the C+G content of the sequence, corresponding to the fraction of (C+G)/(A+T+C+G...). A value of \option{NA} is returned if the function encounters an error that prevents proper calculation of GC percent.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{calcGC-methods}}, \code{\link{calcTm}} }
\examples{
#demo(pipeline,package="HELP")

calcGC("AAAACGCG")
calcGC(sequence="cXgXcXgXcXgX",allow="X")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith }
\keyword{ utilities }
