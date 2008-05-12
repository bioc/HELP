##DEFINE GENERIC FUNCTION calcPrototype()
setGeneric("calcPrototype", 
	function(x, ...) {
		standardGeneric("calcPrototype")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("calcPrototype", 
	signature=c("missing"), 
	function(x, ...) {
		stop("argument 'x' is missing with no defaults")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("calcPrototype", 
	signature=c("ExpressionSet"), 
	function (x, element="exprs", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		callGeneric(assayDataElement(x, element), ...)
	}
)

##DEFINE MAIN calcPrototype() METHOD TO HANDLE CLASS: "matrix"
setMethod("calcPrototype", 
	signature=c("vector"), 
	function(x, samples=NULL, center=TRUE, trim=0.1, verbose=FALSE, ...) {
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Calculating prototype from", length(samples), "samples ...")
		}
		x <- getSamples(as.matrix(x), samples, ...)
		if (dim(x)[1] < 2) {
			stop("'x' must contain at least 2 rows")
		}
		if (dim(x)[2] < 2) {
			stop("need to use at least 2 samples at a time")
		}
		if (center) {
			if (verbose) {
				cat("\n\tCentering data ... ")
			}
			centers <- apply(x, 2, mean, trim=trim, ...)
			x <- t(t(x) - centers)
			if (verbose) {
				cat("FINISHED")
			}
		}
		if (verbose) {
			cat("\n\tApplying means by row ... ")
		}
		x <- apply(x, 1, mean, trim=trim, ...)
		if (center) {
			x <- x + mean(centers, trim=trim, ...)
		}
		if (verbose) {
			cat("FINISHED (", (proc.time()["elapsed"]-start), "s elapsed)\n", sep="")
		}
		return(x)
	}
)

