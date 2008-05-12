\name{readPairs}
\alias{readPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read Nimblegen .pair files }
\description{
  Function to extract data from corresponding files in the Nimblegen .pair format.
}
\usage{
readPairs(x, y, z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the name of the file containing data from signal channel 1. Each line of the file is interpreted as a single data point. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}()}. Tilde-expansion is performed where supported. 
Alternatively, \code{x} can be a readable connection (which will be opened for reading if necessary, and if so closed at the end of the function call).  
\code{x} can also be a complete URL. }
  \item{y}{ the name of the file containing data from signal channel 1. Each line of the file is interpreted as a single data point. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}()}. Tilde-expansion is performed where supported. 
Alternatively, \code{y} can be a readable connection (which will be opened for reading if necessary, and if so closed at the end of the function call).  
\code{y} can also be a complete URL. }
  \item{z}{ object in which to store pair information from files. Can be an \code{ExpressionSet}, in which case pair data will be stored in \code{featureData}. }
  \item{\dots}{ Arguments to be passed to methods (see \code{\link{readPairs-methods}}):
	\item{\code{name}}{ character input specifying a sample name to assign to the data from specified pair files. If \option{NULL} (default), a value will be extracted from the filename specified in \code{x}.}
	\item{\code{element.x}}{ which element of \code{AssayData} (default is \option{"exprs"}) in which to store signal channel 1 data.}
	\item{\code{element.y}}{ which element of \code{AssayData} (default is \option{"exprs2"}) in which to store signal channel 2 data.}
	\item{\code{match.probes}}{ logical specifying whether to match data from pair files by \option{"PROBE\_ID"} with any pre-existing data. if \option{TRUE} (default), order of values will be rearranged so long as there are already \option{"PROBE\_ID"}s specified in \code{featureData}.}
	\item{\code{path}}{ a character vector containing a single full path name to which filenames will be appended. If \option{NULL}, filenames (\code{x} and \code{y}) are treated as is. }
	\item{\code{comment.char}}{ character: a character vector of length one containing a single character or an empty string (default is \option{"\#"}). Use \option{""} to turn off the interpretation of comments altogether. }
	\item{\code{sep}}{ the field separator character (default is \option{"\t"}). Values on each line of the file are separated by this character. If \option{sep = ""} the separator is "white space", that is one or more spaces, tabs, newlines or carriage returns. }
	\item{\code{quote}}{ the set of quoting characters (default is \option{"\""}). To disable quoting altogether, use \code{quote = ""}. See \code{\link{scan}} for the behaviour on quotes embedded in quotes. Quoting is only considered for columns read as character, which is all of them unless \option{colClasses} is specified. }
	\item{\code{eSet}}{ \code{ExpressionSet} input (default is \code{new("ExpressionSet")}) in which to store pair information in \code{assayData} }
	\item{\code{verbose}}{ logical; if \option{TRUE} (default) progress will be output to screen. If \option{FALSE}, no output will be displayed.}
	\item{\code{\dots}}{ other arguments to be passed to \code{read.table}. See \code{\link{read.table}}. }
  }
}
\value{
  Returns an \code{ExpressionSet} filled with \code{assayData} containing matrices of data from both signal channels.
   
  and \code{featureData} containing the following \code{featureColumns}:
  \item{\code{SEQ_ID}}{a vector of characters with container IDs, linking each probe to a parent identifier}
  \item{\code{PROBE_ID}}{a vector of characters containing unique ID information for each probe}

  and \code{phenoData} containing the following \code{sampleColumns}:
  \item{\code{CHIPS}}{a vector of characters with .pair file locations for signal channel 1 data }
  \item{\code{CHIPS2}}{a vector of characters with .pair file locations for signal channel 2 data }
  
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{readPairs-methods}}, \code{\link{read.table}} }
\examples{
#demo(pipeline,package="HELP")

x <- 1:500
y <- rev(x)
data <- sample(8000:10000/1000,size=500)
seqids <- sample(1:50,size=500,replace=TRUE)
cat("#COMMENT\nSEQ_ID\tPROBE_ID\tX\tY\tPM\n",file="./read.pair.test.1")
table.1 <- cbind(seqids,y,x,x,data)
write.table(table.1,file="./read.pair.test.1",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
cat("#COMMENT\nSEQ_ID\tPROBE_ID\tX\tY\tPM\n",file="./read.pair.test.2")
table.2 <- cbind(seqids,y,x,x,rev(data))
write.table(table.2,file="./read.pair.test.2",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
x <- readPairs("./read.pair.test.1","./read.pair.test.2")
c(seqids[1],y[1],data[1],rev(data)[1])
pData(featureData(x))$"SEQ_ID"[1]
pData(featureData(x))$"PROBE_ID"[1]
assayDataElement(x, "exprs")[1]
assayDataElement(x, "exprs2")[1]

#rm(table.1,table.2,x,y,data,seqids)
#file.remove("./read.pair.test.1")
#file.remove("./read.pair.test.2")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
