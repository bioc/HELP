##DEFINE GENERIC FUNCTION readPairs()
setGeneric("readPairs", 
	function(x, y, z, ...) {
		standardGeneric("readPairs")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("readPairs", 
	signature=c("missing", "missing", "missing"), 
	function(x, y, z, ...) {
		stop("argument 'x' is missing with no default")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("readPairs", 
	signature=c("vector", "missing", "missing"), 
	function(x, y, z, ...) {
		if (length(x) > 2) {
			stop("argument 'x' must specify only one or two files")
		}
		if (length(x) <= 1) {
			stop("argument 'y' is missing with no default")
		}
		callGeneric(x[1], x[2], new("ExpressionSet"), ...)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("readPairs", 
	signature=c("vector", "vector", "missing"), 
	function(x, y, z, ...) {
		if (length(x) != 1) {
			stop("argument 'x' must specify one file")
		}
		if (length(y) != 1) {
			stop("argument 'y' must specify one file")
		}
		callGeneric(as.character(x), as.character(y), new("ExpressionSet"), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "character"
setMethod("readPairs", 
	signature=c("character", "character", "ExpressionSet"), 
	function(x, y, z, name=NULL, element.x="exprs", element.y="exprs2", match.probes=TRUE, path=NULL, comment.char="#", sep="\t", quote="\"", verbose=FALSE, ...) {
		if (is.null(name)) {
			name <- sub(paste("[_.]*([0-9][0-9][0-9]|[Cc][Yy][0-9])[_.]*", "[Pp][Aa][Ii][Rr].*", sep=""), "", x)
			name <- sub(".*/", "", name)
		}
		if (!is.null(path)) {
			x <- file.path(path, x)
			y <- file.path(path, y)	
		}
		if (!file.exists(x)) {
			stop("PAIR file does not exist: ", x)
		}
		if (!file.exists(y)) {
			stop("PAIR file does not exist: ", y)
		}
 		if (!validObject(z)) {
			stop("'z' not a valid ExpressionSet object")
		}
 		if (element.x == element.y) {
			stop("'element.x' and 'element.y' must be unique")
		}
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Reading PAIR file:", x, "... ")
		}
		data.x <- read.table(x, as.is=TRUE, header=TRUE, comment.char=comment.char, sep=sep, quote=quote, ...)
		headers.x <- colnames(data.x)
		if (any(is.na(match(c("SEQ_ID", "PROBE_ID", "PM", "X", "Y"), headers.x)))) {
  			stop("one or more required headers (SEQ_ID, PROBE_ID, PM, X, Y) missing")
  		}
 		if (!is.numeric(data.x$"PM")) {
  			stop("'PM' values not numeric in ", x)
    	}
    	if (any(data.x$"PM" <= 0 | data.x$"PM" > 7E4)) {
    		stop("'PM' values out of expected range (0<PM<7E4) in ", x)
    	}
		if (verbose) {
			end <- proc.time()["elapsed"]
			cat("FINISHED (", (end-start), "s elapsed)\n")
			start <- end
			cat("Reading PAIR file:", y, "... ")
		}
		data.y <- read.table(y, as.is=TRUE, header=TRUE, comment.char=comment.char, sep=sep, quote=quote, ...)
		if (verbose) {
			end <- proc.time()["elapsed"]
			cat("FINISHED (", (end-start), "s elapsed)\n")
			start <- end
			cat("Extracting features ... ")
		}
		headers.y <- colnames(data.y)
		if (any(is.na(match(c("SEQ_ID", "PROBE_ID", "PM", "X", "Y"), headers.y)))) {
  			stop("one or more required headers (SEQ_ID, PROBE_ID, PM, X, Y) missing")
  		}
 		if (!is.numeric(data.y$"PM")) {
  			stop("'PM' values not numeric in ", y)
    	}
    	if (any(data.y$"PM" <= 0 | data.y$"PM" > 7E4)) {
    		stop("'PM' values out of expected range (0<PM<7E4) in ", y)
    	}
    	if (dim(data.x)[1] != dim(data.y)[1]) {
    		stop("PAIR file lengths do not match")
    	}
    	data.x$"PROBE_ID" <- paste(data.x$"PROBE_ID", data.x$"X", data.x$"Y", sep="_")
    	data.y$"PROBE_ID" <- paste(data.y$"PROBE_ID", data.y$"X", data.y$"Y", sep="_")
    	if (any(data.x$"PROBE_ID" != data.y$"PROBE_ID")) {
    		data.y <- data.y[match(data.x$"PROBE_ID", data.y$"PROBE_ID"), ]
    		if (any(is.na(data.y$"PROBE_ID"))) {
    			stop("PAIR file probes do not match")
    		}
    	}
		if (match.probes & ("PROBE_ID" %in% featureNames(z))) {
			data.x <- data.x[match(getFeatures(z, "PROBE_ID"), data.x$"PROBE_ID"), ]
			data.y <- data.y[match(getFeatures(z, "PROBE_ID"), data.y$"PROBE_ID"), ]
		}
		else {
			features <- data.frame(cbind(data.x$"SEQ_ID", data.x$"PROBE_ID"), row.names=data.x$"PROBE_ID")
			names(features) <- c("SEQ_ID", "PROBE_ID")
			dimLabels(featureData(z)) <- c("featureNames", "featureColumns")
			featureData(z) <- combine(featureData(z), new("AnnotatedDataFrame", data=features, dimLabels=c("featureNames", "featureColumns")))
		}
		if (verbose) {
			end <- proc.time()["elapsed"]
			cat("FINISHED (", (end-start), "s elapsed)\n")
			start <- end
			cat("Updating ExpressionSet with PAIR data ... ")
		}
		assaydata.x <- assayDataElement(z, element.x)
		assaydata.y <- assayDataElement(z, element.y)
		data.x <- as.matrix(data.x$"PM")
		data.y <- as.matrix(data.y$"PM")
		colnames(data.x) <- colnames(data.y) <- name
		if (is.null(assaydata.x) | all(dim(assaydata.x) == 0)) {
			assaydata.x <- data.x	
		}
		else {
			assaydata.x <- cbind(assaydata.x, data.x)
		}		
		if (is.null(assaydata.y) | all(dim(assaydata.y) == 0)) {
			assaydata.y <- data.y	
		}
		else {
			assaydata.y <- cbind(assaydata.y, data.y)
		}
		assayDataElement(z, element.x) <- assaydata.x
		assayDataElement(z, element.y) <- assaydata.y
		name <- c(sampleNames(z), name)
		phenodata <- pData(phenoData(z))
		empty <- matrix(NA, nrow=dim(phenodata)[1], ncol=1)
		if (!("CHIPS" %in% colnames(phenodata))) {
			colnames(empty) <- "CHIPS"
			phenodata <- cbind(phenodata, empty)
		}
		if (!("CHIPS2" %in% colnames(phenodata))) {
			colnames(empty) <- "CHIPS2"
			phenodata <- cbind(phenodata, empty)
		}
		empty <- matrix(NA, nrow=1, ncol=dim(phenodata)[2])
		colnames(empty) <- colnames(phenodata)
		phenodata <- rbind(phenodata, empty)
		phenodata$"CHIPS"[dim(phenodata)[1]] <- x
		phenodata$"CHIPS2"[dim(phenodata)[1]] <- y
		rownames(phenodata) <- name
		phenoData(z) <- new("AnnotatedDataFrame", data=as.data.frame(phenodata), dimLabels=c("sampleNames", "sampleColumns"))
		protocolData(z) <- phenoData(z)		
		if (verbose) {
			end <- proc.time()["elapsed"]
			cat("FINISHED (", (end-start), "s elapsed)\n")
		}
		if (length(featureNames(assayData(z))) != length(featureNames(featureData(z)))) {
			featureNames(assayData(z)) <- featureNames(featureData(z))
		}
		return(z)
	}
)

