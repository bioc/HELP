##DEFINE GENERIC FUNCTION getSamples()
setGeneric("getSamples", 
	function(x, y, ...) {
		standardGeneric("getSamples")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("getSamples", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getSamples", 
	signature=c("ExpressionSet", "missing"), 
	function(x, y, element="exprs", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!(element %in% assayDataElementNames(x))) {
			stop("'element' not a valid element of ExpressionSet")
		}
		return(as.matrix(assayDataElement(x, element)))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getSamples", 
	signature=c("ExpressionSet", "NULL"), 
	function(x, y, element="exprs", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!(element %in% assayDataElementNames(x))) {
			stop("'element' not a valid element of ExpressionSet")
		}
		return(as.matrix(assayDataElement(x, element)))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getSamples", 
	signature=c("vector", "missing"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getSamples", 
	signature=c("vector", "NULL"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("getSamples", 
	signature=c("matrix", "missing"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("getSamples", 
	signature=c("matrix", "NULL"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getSamples", 
	signature=c("ExpressionSet", "vector"), 
	function(x, y, element="exprs", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!(element %in% assayDataElementNames(x))) {
			stop("'element' not a valid element of ExpressionSet")
		}
		callGeneric(assayDataElement(x, element), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("getSamples", 
	signature=c("vector", "vector"), 
	function(x, y, ...) {
		callGeneric(as.matrix(x), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getSamples", 
	signature=c("matrix", "vector"), 
	function(x, y, order=NULL, ...) {
		x <- as.matrix(x)
		y <- as.vector(y)
		if (is.null(colnames(x))) {
			colnames(x) <- 1:dim(x)[2]
		}
		samples <- colnames(x)
		which.samples <- match(y, samples)
		y[which(!is.na(which.samples))] <- which.samples[which(!is.na(which.samples))]
		y <- as.integer(y)
		if (any(is.na(y))) {
			warning("unmatched samples")
		}
		y <- y[which(!is.na(y))]
		y <- y[which(1<=y & y<=length(samples))]
		if (length(y)<1) {
			stop("no matching samples")	
		}
		y <- unique(y)
		if (length(order) == 1) {
			if (order %in% samples) {
				order <- match(order, samples)
			}
			else {
				order <- as.integer(order)	
			}
			order <- order[which(1<=order & order<=length(samples))]
			if ((is.integer(order)) & (length(order) == 1)) {
				return(as.matrix(x[order(x[, order]), y]))
			}
		}
		
		return(as.matrix(x[, y]))
	}
)
