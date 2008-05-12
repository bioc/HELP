\name{readSampleKey}
\alias{readSampleKey}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read sample key }
\description{
  Function to extract sample key data from a file and link chip ID information with aliases if they exist.
}
\usage{
readSampleKey(file = NULL, chips = NULL, comment.char = "#", sep = "\t")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ the name of the file containing sample key information. Each line of the file is interpreted as a single chip-to-sample map. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link{getwd}()}. Tilde-expansion is performed where supported. 
Alternatively, \code{file} can be a readable connection (which will be opened for reading if necessary, and if so closed at the end of the function call).  
\option{file} can also be a complete URL. }
  \item{chips}{ a character vector specifying a specific chip ID lookup in the sample key, for which the function will return the appropriate sample aliases }
  \item{comment.char}{ character: a character vector of length one containing a single character or an empty string. Use \option{""} to turn off the interpretation of comments altogether. }
  \item{sep}{ the field separator character. Values on each line of the file are separated by this character. If \option{sep = ""} the separator is "white space", that is one or more spaces, tabs, newlines or carriage returns. }
}
\value{
  Returns a character vector of sample alias information corresponding to the chips present in the sample key or a subset thereof, specified by the \code{chips} input.
}
\author{ Reid F. Thompson (\email{rthompso@aecom.yu.edu}) }
\seealso{ \code{\link{read.table}} }
\examples{
#demo(pipeline,package="HELP")

cat("#COMMENT\nCHIP_ID\tSAMPLE\n",file="./sample.key.txt")
write.table(cbind(1:10,1001:1010),file="./sample.key.txt",append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
readSampleKey(file="./sample.key.txt")
readSampleKey(file="./sample.key.txt",chips=c(7:10,"NA1","NA2"))

#file.remove("./sample.key.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ IO }
