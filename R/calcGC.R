##DEFINE GENERIC FUNCTION calcGC()
setGeneric("calcGC", 
	function(x, ...) {
		standardGeneric("calcGC")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("calcGC", 
	signature=c("missing"), 
	function(x, ...) {
		return(NULL)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "NULL"
setMethod("calcGC", 
	signature=c("NULL"), 
	function(x, ...) {
		return(NULL)
	}
)



##DEFINE METHOD TO HANDLE CLASS: "character"
setMethod("calcGC", 
	signature=c("character"), 
	function(x, allow="N", ...) {
		x <- toupper(x)
		allow <- toupper(allow)
		seq.lengths <- nchar(x)
		sequence <- rep(NA,length(x))
		chars.all <- strsplit(x,split="")
		for (i in 1:length(chars.all)) {
			chars <- chars.all[[i]]
			seq.length <- seq.lengths[i]
			if (seq.length <= 0) {
				warning("Unable to calculate G+C content: empty sequence")
				next			
			}
			if (any(!(chars %in% c("A", "T", "C", "G", allow)))) {
				warning("Unable to calculate G+C content: unrecognized characters")
				next
			}
			sequence[i] <- length(which(chars == "C" | chars == "G"))/seq.length
		}
		return(sequence)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("calcGC", 
	signature=c("ExpressionSet"), 
	function(x, ...) {
		if (!validObject(x)) {
			stop("'x' not a valid ExpressionSet object")
		}
		if (!"SEQUENCE" %in% varLabels(featureData(x))) {
			stop("featureData must contain 'SEQUENCE' information")
		}
		if ("GC" %in% varLabels(featureData(x))) {
			stop("featureData already contains GC content information")
		}
		sequence <- as.character(pData(featureData(x))[, "SEQUENCE"])
		gc <- callGeneric(sequence, ...)
		gc <- data.frame(GC=gc)
		featureData(x) <- combine(featureData(x),new("AnnotatedDataFrame", data=gc, dimLabels=c("featureNames", "featureColumns")))
		return(x)
	}
)


