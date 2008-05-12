\name{calcTm}
\alias{calcTm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Tm }
\description{
  Calculate melting temperature (Tm) using the nearest-neighbor base-stacking algorithm and the unified thermodynamic parameters.
}
\usage{
calcTm(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ characters containing nucleotide sequences (ex: \option{"ATCGGAA"}) or an object of class \code{ExpressionSet} }
  \item{\dots}{ Additional arguments passed to methods:
	\code{\item{strand1.concentration}}{ numeric value specifying concentration of strand 1 (default is 2e-07)}
	\code{\item{strand2.concentration}}{ numeric value specifying concentration of strand 2 (default is 2e-07)}
	\code{\item{method}}{ character value specifying the Tm algorithm to use (default is \option{"nearest-neighbor"}); currently not supported }
  }
}
\value{
  Returns a numerical value indicating the predicted melting temperature (Tm) of the sequence in degrees Celsius. A value of \option{NA} is returned if the function encounters an error that prevents proper Tm calculation.
}
\references{   Allawi, H.T. and SantaLucia, J., Jr. (1997) Thermodynamics and NMR of internal G.T mismatches in DNA, \emph{Biochemistry}, \strong{36}, 10581-10594. }
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{calcTm-methods}}, \code{\link{base.stacking.thermodynamics}}, \code{\link{calcGC}} }
\examples{
#demo(pipeline,package="HELP")

calcTm("GTGTGGCTACAGGTGGGCCGTGGCGCACCTAAGTGAGGACAGAGAACAAC")
calcTm("GTGTGGCTACAGGTGGGCCGTGGCGCACCTAAGTGAGGACAGAGAACAAC",strand1.concentration=1E-5,strand2.concentration=2E-8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ arith }
\keyword{ utilities }
