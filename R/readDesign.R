##DEFINE GENERIC FUNCTION readDesign()
setGeneric("readDesign", 
	function(x, y, z, ...) {
		standardGeneric("readDesign")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("readDesign", 
	signature=c("missing", "missing", "missing"), 
	function(x, y, z, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("readDesign", 
	signature=c("vector", "missing", "missing"), 
	function(x, y, z, ...) {
		x <- as.vector(x)
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
setMethod("readDesign", 
	signature=c("vector", "vector", "missing"), 
	function(x, y, z, ...) {
		x <- as.vector(x)
		y <- as.vector(y)
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
setMethod("readDesign", 
	signature=c("character", "character", "ExpressionSet"), 
	function(x, y, z, path=NULL, comment.char="#", sep="\t", quote="\"", verbose=FALSE, ...) {
		if (!is.null(path)) {
			x <- file.path(path, x)	
			y <- file.path(path, y)	
		}
		if (!file.exists(x)) {
			stop("NDF file does not exist: ", x)
		}
		if (!file.exists(y)) {
			stop("NGD file does not exist: ", y)
		}
		if (!validObject(z)) {
			stop("'z' not a valid ExpressionSet object")
		}
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Reading NDF file:", x, "... ")	
		}
		ndf <- read.table(x, as.is=TRUE, header=TRUE, comment.char=comment.char, sep=sep, quote=quote, ...)
		headers.ndf <- colnames(ndf)
  		if (any(is.na(match(c("SEQ_ID", "X", "Y", "PROBE_ID", "CONTAINER", "PROBE_SEQUENCE", "PROBE_DESIGN_ID"), headers.ndf)))) {
  			stop("one or more required headers (SEQ_ID, PROBE_ID, X, Y, CONTAINER, PROBE_SEQUENCE, PROBE_DESIGN_ID) missing from ", x)
  		}
    	ndf$"PROBE_ID" <- paste(ndf$"PROBE_ID", ndf$"X", ndf$"Y", sep="_")
		if (verbose) {
			cat("FINISHED\nReading NGD file:", y, "... ")
		}
		ngd <- read.table(y, as.is=TRUE, header=TRUE, comment.char=comment.char, sep=sep, quote=quote, ...)
		headers.ngd <- colnames(ngd)
	  	if (any(is.na(match(c("SEQ_ID", "CHROMOSOME", "START", "STOP"), headers.ngd)))) {
	  		stop("one or more required headers (SEQ_ID, CHROMOSOME, START, STOP) missing from ", y)
	  	}
		if (verbose) {
			cat("FINISHED\nProcessing designs ... ")
		}
		control.probes  <- c(grep("CONTROLS", ndf$"CONTAINER"), grep("_CODE", ndf$"CONTAINER"))
		random.probes  <- grep("RAND", ndf$"PROBE_ID")
		real.probes <- which(!c(1:length(ndf$"SEQ_ID")) %in% c(control.probes, random.probes))
		if (any(is.na(match(unique(ndf$"SEQ_ID"[real.probes]), ngd$"SEQ_ID")))) {
			stop("one or more SEQ_ID mismatches: missing from NGD")
		}
		if (any(is.na(match(ngd$"SEQ_ID", unique(ndf$"SEQ_ID"))))) {
			stop("one or more SEQ_ID mismatches: missing from NDF")
		}
		ndf$"CONTAINER"[random.probes] <- "RAND"
		rm(random.probes)
		ndf$"CONTAINER"[control.probes] <- "CONTROL"
		ndf$"CONTAINER"[real.probes] <- "DATA"
		real.probes <- which(!c(1:length(ndf$"SEQ_ID")) %in% c(control.probes))
		ndf <- ndf[real.probes, ]
		rm(real.probes, control.probes)
		if (verbose) {
			cat("FINISHED\nFormatting information ... ")
		}
		feature.names <- featureNames(featureData(z))
		if (length(feature.names) > 1) {
			feature.ord <- match(feature.names, ndf$"PROBE_ID")
			ndf <- ndf[feature.ord, ]
		}
		map <- match(ndf$"SEQ_ID", ngd$"SEQ_ID")
		features <- data.frame(cbind(ndf$"SEQ_ID", ndf$"PROBE_ID", ndf$"X", ndf$"Y", ndf$"CONTAINER", ngd$"CHROMOSOME"[map], ngd$"START"[map], ngd$"STOP"[map], as.numeric((ngd$"STOP"-ngd$"START"+1)[map]), ndf$"PROBE_SEQUENCE"), row.names=1:dim(ndf)[1])
		names(features) <- c("SEQ_ID", "PROBE_ID", "X", "Y", "TYPE", "CHR", "START", "STOP", "SIZE", "SEQUENCE")
		rm(ngd)
		if ("DMD" %in% headers.ndf) {
			map <- as.matrix(ndf$"DMD")
			colnames(map) <- "WELL"			
			features <- cbind(features, map)
		}
		rm(ndf, map)
		featureData(z) <- new("AnnotatedDataFrame", data=features, dimLabels=c("featureNames", "featureColumns"))
		if (verbose) {
			cat("FINISHED (", (proc.time()["elapsed"]-start), "s elapsed)\n")
		}
		return(z)
	}
)

