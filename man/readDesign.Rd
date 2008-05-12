\name{readDesign}
\alias{readDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read NimbleGen design files }
\description{
  Function to extract array design information from corresponding files in the Nimblegen .ndf and .ngd formats.
}
\usage{
readDesign(x, y, z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ path to the Nimblegen design file (.ndf). Each line of the file is interpreted as a single spot on the array design. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}()}. Tilde-expansion is performed where supported. 
Alternatively, \code{x} can be a readable connection (which will be opened for reading if necessary, and if so closed at the end of the function call).  
\option{file} can also be a complete URL. }
  \item{y}{ path to the Nimblegen gene descriptions file (.ngd). Each line of the file is interpreted as a single locus. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}()}. Tilde-expansion is performed where supported. 
Alternatively, \code{y} can be a readable connection (which will be opened for reading if necessary, and if so closed at the end of the function call).  
\option{file} can also be a complete URL. }
  \item{z}{ object in which to store design information from files. Can be an \code{ExpressionSet}, in which case design information will be stored in \code{featureData}.  }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{readDesign-methods}}):
	\item{\code{path}}{ a character vector containing a single full path name to which filenames will be appended. If \option{NULL}, filenames (\code{x} and \code{y}) are treated as is. }
	\item{\code{comment.char}}{ character: a character vector of length one containing a single character or an empty string (default is \option{"\#"}). Use \option{""} to turn off the interpretation of comments altogether. }
	\item{\code{sep}}{ the field separator character (default is \option{"\t"}). Values on each line of the file are separated by this character. If \option{sep = ""} the separator is "white space", that is one or more spaces, tabs, newlines or carriage returns. }
	\item{\code{quote}}{ the set of quoting characters (default is \option{"\""}). To disable quoting altogether, use \code{quote = ""}. See \code{\link{scan}} for the behavior on quotes embedded in quotes. Quoting is only considered for columns read as character, which is all of them unless \option{colClasses} is specified. }
	\item{\code{eSet}}{ \code{ExpressionSet} input (default is \code{new("ExpressionSet")}) in which to store design information in \code{featureData} }
	\item{\code{\dots}}{ other arguments to be passed to \code{read.table}. See \code{\link{read.table}}. }
  }
}
\value{
  Returns an \code{ExpressionSet} filled with \code{featureData} containing the following \code{featureColumns}:
  \item{\code{SEQ_ID}}{a vector of characters with container IDs, linking each probe to a parent identifier}
  \item{\code{PROBE_ID}}{a vector of characters containing unique ID information for each probe}
  \item{\code{X}}{vector of numerical data determining x-coordinates of probe location on chip}
  \item{\code{Y}}{vector of numerical data determining y-coordinates of probe location on chip}
  \item{\code{TYPE}}{a vector of characters defining the type of probe, e.g. random background signals (\option{"RAND"}) or usable data (\option{"DATA"}).}
  \item{\code{CHR}}{a matrix of characters containing unique ID and chromosomal positions for each container}
  \item{\code{START}}{a matrix of characters containing unique ID and chromosomal positions for each container}
  \item{\code{STOP}}{a matrix of characters containing unique ID and chromosomal positions for each container}
  \item{\code{SIZE}}{a matrix of characters containing unique ID and chromosomal positions for each container}
  \item{\code{SEQUENCE}}{a vector of characters containing sequence information for each probe}
  \item{\code{WELL}}{a vector of characters containing multiplex well location for each probe (if present in design files)}

}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{readDesign-methods}}, \code{\link{read.table}} }
\examples{
#demo(pipeline, package="HELP")

chr <- rep("chr1", 500)
start <- (1:500)*200
stop <- start+199
x <- 1:500
seqids <- sample(1:50, size=500, replace=TRUE)
cat("#COMMENT\nSEQ_ID\tCHROMOSOME\tSTART\tSTOP\n", file="./read.design.test.ngd")
table.ngd <- cbind(seqids, chr, start, stop)
write.table(table.ngd, file="./read.design.test.ngd", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
cat("#COMMENT\nSEQ_ID\tX\tY\tPROBE_ID\tCONTAINER\tPROBE_SEQUENCE\tPROBE_DESIGN_ID\n", file="./read.design.test.ndf")
sequence <- rep("NNNNNNNN", 500)
table.ndf <- cbind(seqids, x, x, x, x, sequence, x)
write.table(table.ndf, file="./read.design.test.ndf", append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
x <- readDesign("./read.design.test.ndf", "./read.design.test.ngd")
seqids[1:10]
pData(featureData(x))$"SEQ_ID"[1:10]

#rm(table.ngd, table.ndf, chr, start, stop, x, seqids, sequence)
#file.remove("./read.design.test.ngd")
#file.remove("./read.design.test.ndf")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
