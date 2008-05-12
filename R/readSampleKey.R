readSampleKey <- function(file=NULL,chips=NULL,comment.char="#",sep="\t") {
	if (!file.exists(file)) {
		warning("Cannot open sample key ",file)
		return(chips)
	}
	sample.key <- read.table(file,header=TRUE,as.is=TRUE,sep=sep,row.names=NULL,comment.char=comment.char)
	headers <- colnames(sample.key)
	if (any(is.na(match(c("CHIP_ID","SAMPLE"),headers)))) {
		warning("one or more required headers (CHIP_ID, SAMPLE) missing from ",file)
		return(chips)
	}
	if (is.null(chips)) {
		chips <- as.character(sample.key$"SAMPLE")
	}
	chip.samples <- as.character(sample.key$"SAMPLE"[match(chips,as.character(sample.key$"CHIP_ID"))])
	not.na <- which(!is.na(chip.samples))
	chips[not.na] <- chip.samples[not.na]
	return(chips)
}
