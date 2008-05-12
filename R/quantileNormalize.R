##DEFINE GENERIC FUNCTION quantileNormalize()
setGeneric("quantileNormalize", 
	function(x, y, ...) {
		standardGeneric("quantileNormalize")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("quantileNormalize", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("quantileNormalize", 
	signature=c("matrix", "missing"), 
	function(x, y, ...) {
		if (dim(x)[2] != 2) {
			stop("argument 'y' is missing with no default")
		}
		callGeneric(x[, 1], y=x[, 2], ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("quantileNormalize", 
	signature=c("vector", "missing"), 
	function(x, y, ...) {
		stop("argument 'y' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("quantileNormalize", 
	signature=c("ExpressionSet","missing"), 
	function (x, y, element="exprs", sample=1, feature="SIZE", ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), getFeatures(x, feature), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("quantileNormalize", 
	signature=c("ExpressionSet","vector"), 
	function (x, y, element="exprs", sample=1, ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("quantileNormalize", 
	signature=c("vector", "ExpressionSet"), 
	function (x, y, feature="SIZE", ...) {
		if (!validObject(y)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(x, getFeatures(y, feature), ...)
	}
)
	
##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("quantileNormalize", 
	signature=c("ExpressionSet","ExpressionSet"), 
	function (x, y, element="exprs", sample=1, feature=1, ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!validObject(y)) {
			stop("argument 'y' not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), getFeatures(y, feature), ...)
	}
)


##DEFINE MAIN quantileNormalize() FUNCTION
setMethod("quantileNormalize",
	signature=c("vector","vector"), 
	function (x, y, num.bins=10, num.steps=3, mode=c("continuous", "discrete"), type=7, na.rm=TRUE, ...) {
		mode <- match.arg(mode)
		rows <- rownames(x)
		x <- as.numeric(x)
		y <- as.numeric(y)
		if (length(x) != length(y)) {
			stop("'x' and 'y' lengths do not match")
		}
		if (any(is.na(x))) {
			if (na.rm) {
				not.na <- which(!is.na(x))
				y <- y[not.na]
				x <- x[not.na]
				if (!is.null(rows)) {
					rows <- rows[not.na]
				}		
			}
			else {
				stop("'x' contains missing values")
			}
		}
		if (any(is.na(y))) {
			if (na.rm) {
				not.na <- which(!is.na(y))
				y <- y[not.na]
				x <- x[not.na]		
				if (!is.null(rows)) {
					rows <- rows[not.na]
				}		
			}
			else {
				stop("'y' contains missing values")
			}
		}
		length.x <- length(x)
		discrete.y <- unique(y)
		length.y <- length(discrete.y)	
		num.steps <- max(as.integer(num.steps),1)
		if (num.steps > length.y) {
			stop("'num.steps' cannot exceed number of data points")
		}
		num.bins <- max(as.integer(num.bins),1)
		if (num.bins > length.y) {
			stop("'num.bins' cannot exceed number of data points")
		}
		num.total <- num.steps*(num.bins-1)+1
		dim.quantiles <- 2*length(y)/num.bins
		bounds <- matrix(nrow=num.bins+1, ncol=num.steps)
		if (mode == "continuous") {
			for (step in 1:num.steps) {
				bounds[1:num.bins, step] <- quantile(y, ((0:(num.bins-1))+(step-1)/num.steps)/num.bins)
				bounds[num.bins+1, step] <- max(y)
			}
		}
		if (mode == "discrete") {
			for (step in 1:num.steps) {
				bounds[1:num.bins, step] <- quantile(discrete.y, ((0:(num.bins-1))+(step-1)/num.steps)/num.bins)
				bounds[num.bins+1, step] <- max(discrete.y)
			}
		}
		rm(discrete.y)
		IDs <- quantiles <- matrix(nrow=dim.quantiles, ncol=num.total)
		for (step in 1:num.steps) {
			for ( i in 1:num.bins) {
				bin <- (i-1)*num.steps+step
				if (bin > num.total) {
					next
				}
				if (bounds[i, step] > bounds[i+1, step]) {
					bounds[i+1, step] <- bounds[i, step]
				}
				bins <- which(y >= bounds[i, step] & y <= bounds[i+1, step])
				if (length(bins) < 2) {
					warning("Skipping bin ", bin, " (step ", step, ", bin ", i, "): need at least 2 data points")
					next
				}
  				quantiles[, bin] <- quantile(x[bins], 0:(dim.quantiles-1)/(dim.quantiles-1), na.rm=na.rm, type=type, ...)
  				if (length(bins) > dim(IDs)[1]) {
					warning("Excluding data for bin ", bin, " (step ", step, ", bin ", i, "): too many points, analyzing a random subset")
  					bins <- sample(bins, size=dim(IDs)[1])
  				}
   				IDs[1:length(bins), bin] <- bins[order(x[bins])]
			}
		}
		average <- apply(quantiles, 1, mean, na.rm=na.rm)
		quantiles <- matrix(nrow=length(x), ncol=2, data=0)
		for (step in 1:num.steps) {
			for ( i in 1:num.bins) {
				bin <- (i-1)*num.steps+step
				if (bin > num.total) {
					next
				}
				bin.IDs <- IDs[which(!is.na(IDs[, bin])), bin]
				quantiles[bin.IDs, 1] <- quantiles[bin.IDs, 1]+quantile(average, 0:(length(bin.IDs)-1)/(length(bin.IDs)-1), na.rm=na.rm, type=type, ...)
				quantiles[bin.IDs, 2] <- quantiles[bin.IDs, 2]+1
			}
		}
		bins <- which(quantiles[, 2]>0)
		x[bins] <- quantiles[bins, 1]/quantiles[bins, 2]
		if (!is.null(rows)) {
			x <- as.matrix(x)
			rownames(x) <- rows
		}
		return(x)
	}
)
