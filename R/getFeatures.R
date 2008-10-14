##DEFINE GENERIC FUNCTION getFeatures()
setGeneric("getFeatures", 
	function(x, y, ...) {
		standardGeneric("getFeatures")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("getFeatures", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("arguments 'x' and 'y' are missing with no defaults")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getFeatures", 
	signature=c("ExpressionSet", "missing"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		x <- featureData(x)
		if (length(varLabels(x)) < 1) {
			warning("empty featureData in ExpressionSet")
			return(matrix(nrow=0, ncol=0))	
		}
		return(as.matrix(pData(x)))
	}
)


##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getFeatures", 
	signature=c("ExpressionSet", "NULL"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		x <- featureData(x)
		if (length(varLabels(x)) < 1) {
			warning("empty featureData in ExpressionSet")
			return(matrix(nrow=0, ncol=0))	
		}
		return(as.matrix(pData(x)))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "AnnotatedDataFrame"
setMethod("getFeatures", 
	signature=c("AnnotatedDataFrame", "missing"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid AnnotatedDataFrame object")
		}
		return(as.matrix(pData(x)))
	}
)


##DEFINE METHOD TO HANDLE CLASS: "AnnotatedDataFrame"
setMethod("getFeatures", 
	signature=c("AnnotatedDataFrame", "NULL"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid AnnotatedDataFrame object")
		}
		return(as.matrix(pData(x)))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getFeatures", 
	signature=c("vector", "missing"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getFeatures", 
	signature=c("vector", "NULL"), 
	function(x, y, ...) {
		return(as.matrix(x))
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("getFeatures", 
	signature=c("ExpressionSet", "vector"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		x <- featureData(x)
		if (length(varLabels(x)) < 1) {
			warning("empty featureData in ExpressionSet")
			return(matrix(nrow=0, ncol=0))	
		}
		callGeneric(as.matrix(pData(x)), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "AnnotatedDataFrame"
setMethod("getFeatures", 
	signature=c("AnnotatedDataFrame", "vector"), 
	function(x, y, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid AnnotatedDataFrame object")
		}
		callGeneric(as.matrix(pData(x)), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("getFeatures", 
	signature=c("vector", "vector"), 
	function(x, y, ...) {
		callGeneric(as.matrix(x), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("getFeatures", 
	signature=c("matrix", "vector"), 
	function(x, y, order=NULL, ...) {
		x <- as.matrix(x)
		y <- as.vector(y)
		if (is.null(colnames(x))) {
			colnames(x) <- 1:dim(x)[2]
		}
		features <- colnames(x)
		y <- fuzzyMatches(y, colnames(x), strict=FALSE, keep=FALSE, alias=FALSE)
#		which.features <- match(y,features)
		order <- fuzzyMatches(order, colnames(x), strict=FALSE, keep=FALSE, alias=FALSE)
#		y[which(!is.na(which.features))] <- which.features[which(!is.na(which.features))]
		y <- as.numeric(y)
		if (length(y) < 1) {
			stop("no matching features")	
		}
		if (!is.null(order)) {
			return(as.matrix(x[order(x[, order]), y]))
		}
		return(as.matrix(x[, y]))
	}
)
