##DEFINE GENERIC FUNCTION combineData()
setGeneric("combineData", 
	function(x, y, w, ...) {
		standardGeneric("combineData")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("combineData", 
	signature=c("missing", "missing", "missing"), 
	function(x, y, w, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("combineData", 
	signature=c("vector", "missing", "missing"), 
	function(x, y, w, ...) {
		x <- as.matrix(x)
		callGeneric(x, rep(1, dim(x)[1]), matrix(data=1, nrow=dim(x)[1], ncol=dim(x)[2]), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("combineData", 
	signature=c("vector", "missing", "vector"), 
	function(x, y, w, ...) {
		x <- as.matrix(x)
		callGeneric(x, rep(1, dim(x)[1]), w, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("combineData", 
	signature=c("vector", "vector", "missing"), 
	function(x, y, w, ...) {
		x <- as.matrix(x)
		if (length(y) != dim(x)[1]) {
			callGeneric(x, rep(1, dim(x)[1]), y, ...)
		}
		else {
			callGeneric(x, y, matrix(data=1, nrow=dim(x)[1], ncol=dim(x)[2]), ...)	
		}
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("combineData", 
	signature=c("ExpressionSet", "missing", "missing"), 
	function(x, y, w, element="exprs", feature.group=NULL, element.weight=NULL, feature.weight=NULL, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!is.null(feature.group)) {
			y <- getFeatures(x, feature.group)
		}
		else {
			y <- rep(1, dims(x)[1])	
		}
		if (!is.null(element.weight)) {
			w <- assayDataElement(x, element.weight)
		}
		else {
			if (!is.null(feature.weight)) {
				w <- getFeatures(x, feature.weight)
			}
			else {
				w <- matrix(1, nrow=dims(x)[1], ncol=dims(x)[2])
			}
		}
		callGeneric(assayDataElement(x, element), y, w, ...)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("combineData", 
	signature=c("ExpressionSet", "vector", "missing"), 
	function(x, y, w, element="exprs", element.weight=NULL, feature.weight=NULL, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!is.null(element.weight)) {
			w <- assayDataElement(x, element.weight)
		}
		else {
			if (!is.null(feature.weight)) {
				w <- getFeatures(x, feature.weight)
			}
			else {
				w <- matrix(1, nrow=dims(x)[1], ncol=dims(x)[2])
			}
		}
		callGeneric(assayDataElement(x, element), y, w, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("combineData", 
	signature=c("vector", "vector", "vector"), 
	function(x, y, w, samples=NULL, trim=0, na.rm=FALSE, verbose=FALSE, ...) {
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Combining data (")	
		}
		x <- as.matrix(x)
		y <- as.vector(y)
		w <- as.matrix(w)
		if (is.null(colnames(w))) {
			colnames(w) <- colnames(x)	
		}
		x <- getSamples(x, samples)
		w <- getSamples(w, samples)
		num.obs <- dim(x)[1]
		N <- dim(x)[2]
		if (verbose) {
			cat(N, "samples) ... ")	
		}
		if (num.obs!=length(y)) {
			stop("'x' and 'y' lengths do not match")
		}
		if (any(dim(x) != dim(w))) {
			if (length(w) == num.obs) {
				w <- matrix(data=as.vector(w), nrow=num.obs, ncol=N)
			}
			else {
				stop("'x' and 'w' dimensions do not match")
			}
		}
		if (any(is.na(x))) {
			if (na.rm) {
				w[which(is.na(x))] <- 0
			}
			else {
				stop("'x' contains missing values")	
			}
		}
		trim <- min(0.5, max(0, as.numeric(trim)))
		groups <- unique(y)
		x.combined <- matrix(NA, nrow=length(groups), ncol=N)
		rownames(x.combined) <- groups
		colnames(x.combined) <- colnames(x)
		if (is.null(colnames(x.combined))) {
			colnames(x.combined) <- 1:N
		}
		w[which(w<0)] <- 0
		x.weighted <- x*w
		for (i in 1:length(groups)) {
			y.i <- which(y == groups[i])
			num <- length(y.i)
			if (num == 1) {
				x.combined[i, ] <- x[y.i, ]
				next
			}
			to.trim <- floor(trim*(num-1))
			use.trim <- (to.trim > 0)
			to.trim <- c(1:to.trim, (num+1-to.trim):num)
			for (j in 1:N) {
				if (use.trim) {
					x.order <- y.i[order(x[y.i, j])[to.trim]]
					x.weighted[x.order, j] <- w[x.order, j] <- 0
				}
				x.combined[i, j] <- sum(x.weighted[y.i, j], na.rm=na.rm)/sum(w[y.i, j], na.rm=na.rm)
			}
		}
		if (verbose) {
			cat("FINISHED (", (proc.time()["elapsed"]-start), "s elapsed)\n")
		}
		return(x.combined)
	}
)
