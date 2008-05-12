##DEFINE GENERIC FUNCTION plotChip()
setGeneric("plotChip", 
	function(x, y, z, ...) {
		standardGeneric("plotChip")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("plotChip", 
	signature=c("missing", "missing", "missing"), 
	function(x, y, z, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("plotChip", 
	signature=c("matrix", "missing", "missing"), 
	function(x, y, z, ...) {
		dim.x <- dim(x)[2]
		dim.y <- dim(x)[1]
		x <- x[dim.y:1, ]
		callGeneric(rep(1:dim.x, each=dim.y), rep(1:dim.y, dim.x), x, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotChip", 
	signature=c("ExpressionSet", "missing", "missing"), 
	function(x, y, z, element="exprs", sample=1, feature.x="X", feature.y="Y", ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getFeatures(x, feature.x), getFeatures(x, feature.y), getSamples(x, sample, element=element), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotChip", 
	signature=c("vector", "vector", "ExpressionSet"), 
	function(x, y, z, element="exprs", sample=1, ...) {
		if (!validObject(z)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(x, y, getSamples(z, sample, element=element), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotChip", 
	signature=c("ExpressionSet", "vector", "missing"), 
	function(x, y, z, feature.x="X", feature.y="Y", ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getFeatures(x, feature.x), getFeatures(x, feature.y), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotChip", 
	signature=c("ExpressionSet", "ExpressionSet", "missing"), 
	function(x, y, z, element="exprs", sample=1, feature.x="X", feature.y="Y", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!validObject(y)) {
			stop("argument 'y' not a valid ExpressionSet object")
		}
		callGeneric(getFeatures(x, feature.x), getFeatures(x, feature.y), getSamples(y, sample, element=element), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("plotChip", 
	signature=c("vector", "vector", "vector"), 
	function(x, y, z, na.rm=FALSE, main=NULL, xlab="", ylab="", colors=rev(rainbow(n=20, start=0, end=1/3)), range=c(NA, NA), nrows=NULL, ncols=NULL, ...) {
		x <- as.numeric(x)
		y <- as.numeric(y)
		z <- as.numeric(z)
		range <- as.vector(range)
		colors <- as.vector(colors)
		num.colors <- length(colors)
		if (num.colors < 2) {
			stop("'colors' must have at least 2 values")
		}
		num.obs <- length(x)
		if (num.obs != length(y)) {
			stop("'x' and 'y' lengths do not match")
		}
		if (num.obs != length(z)) {
			stop("'x' and 'z' lengths do not match")
		}
		if (any(is.na(x) | is.na(y) | is.na(z))) {
			if (na.rm) {
				not.na <- which(!is.na(x) & !is.na(y) & !is.na(z))
				x <- x[not.na]
				y <- y[not.na]
				z <- z[not.na]
				rm(not.na)
			}
			else {
				stop("missing values where data expected")
			}	
		}
		if ((length(range) != 2) | (range[1] >= range[2]) | (any(is.na(range)))) {
			range <- c(0, ceiling(max(abs(z))))
			range[1] <- -range[2]
		}
		bounds <- quantile(range, 1:num.colors/num.colors)
		num.cols <- length(unique(x))
		num.rows <- length(unique(y))
		if (is.numeric(nrows) & is.numeric(ncols)) {
			ncols <- ceiling(max(1, min(num.cols,ncols)))
			nrows <- ceiling(max(1, min(num.rows,nrows)))
			x.bounds <- quantile(unique(x), 0:(ncols-1)/ncols)
			x.bounds[ncols+1] <- Inf
			y.bounds <- quantile(unique(y), 0:(nrows-1)/nrows)
			y.bounds[nrows+1] <- Inf
			for (i in 1:ncols) {
				for (j in 1:nrows) {
					block <- which(x >= x.bounds[i] & x < x.bounds[i+1] & y >= y.bounds[j] & y < y.bounds[j+1])
					z[block] <- mean(z[block], na.rm=na.rm)
				}
			}
		}
		z.order <- order(z)
		x <- x[z.order]
		y <- y[z.order]
		z <- z[z.order]
		rm(z.order)
		plot(0, xlim=c(0, max(x)), ylim=c(0, max(y)), type="n", xlab=xlab, ylab=ylab, main=main, ...)
		left <- 1
		for (i in 1:num.colors) {
			right <- which.min(abs(z-bounds[i]))
			if (z[right] > bounds[i]) {
				right <- right-1	
			}
			if (left > right) {
				next	
			}
			rect(x[left:right]-1, y[left:right]-1, x[left:right], y[left:right], border=NA, col=colors[i])
			left <- right+1
		}
		if (num.obs > left) {
			rect(x[left:num.obs]-1, y[left:num.obs]-1, x[left:num.obs], y[left:num.obs], border=NA, col=colors[i])
		}
	}
)

