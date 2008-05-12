##DEFINE GENERIC FUNCTION plotBins()
setGeneric("plotBins", 
	function(x, y, ...) {
		standardGeneric("plotBins")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("plotBins", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("plotBins", 
	signature=c("matrix", "missing"), 
	function(x, y, ...) {
		if (dim(x)[2]!=2) {
			stop("argument 'y' is missing with no default")
		}
		callGeneric(x[, 1], y=x[, 2], ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("plotBins", 
	signature=c("vector", "missing"), 
	function(x, y, ...) {
		stop("argument 'y' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotBins", 
	signature=c("ExpressionSet","missing"), 
	function (x, y, element="exprs", sample=1, feature="SIZE", ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), getFeatures(x, feature), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotBins", 
	signature=c("ExpressionSet","vector"), 
	function (x, y, element="exprs", sample=1, ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotBins", 
	signature=c("vector", "ExpressionSet"), 
	function (x, y, feature=1, ...) {
		if (!validObject(y)) {
			stop("not a valid ExpressionSet object")
		}
		callGeneric(x, getFeatures(y, feature), ...)
	}
)
	
##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotBins", 
	signature=c("ExpressionSet","ExpressionSet"), 
	function (x, y, element="exprs", sample=1, feature="SIZE", ...) {
		if (!validObject(x)) {
			stop("argument 'x' not a valid ExpressionSet object")
		}
		if (!validObject(y)) {
			stop("argument 'y' not a valid ExpressionSet object")
		}
		callGeneric(getSamples(x, sample, element=element), getFeatures(y, feature), ...)
	}
)

##DEFINE MAIN plotBins() FUNCTION
setMethod("plotBins", 
	signature=c("vector","vector"), 
	function (x, y, num.bins=10, num.steps=3, mode=c("continuous", "discrete"), show.avg=FALSE, main=NULL, xlab="", ylab="", na.rm=TRUE, ...) {
		mode <- match.arg(mode)
		x <- as.vector(x)
		y <- as.vector(y)
		if (length(x)!=length(y)) {
			stop("'x' and 'y' lengths do not match")
		}
		if (any(is.na(x))) {
			if (na.rm) {
				not.na <- which(!is.na(x))
				y <- y[not.na]
				x <- x[not.na]		
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
			}
			else {
				stop("'y' contains missing values")
			}
		}
		length.x <- length(x)
		discrete.y <- unique(y)
		length.y <- length(discrete.y)	
		num.steps <- max(as.integer(num.steps), 1)
		if (num.steps > length.y) {
			stop("'num.steps' cannot exceed number of data points")
		}
		num.bins <- max(as.integer(num.bins),1)
		if (num.bins > length.y) {
			stop("'num.bins' cannot exceed number of data points")
		}
		num.total <- num.steps*(num.bins-1)+1
		colors <- rainbow(n=num.total, start=0, end=2/3)
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
		layout(cbind(2, 1), widths=c(6+ceiling(log(num.steps+1)), ceiling(log(num.steps+1))), heights=rep(1, 2))
		plot(c(0, num.steps), c(min(y), max(y)), type="p", cex=0, main="", xlab="", ylab=paste(num.total, " sliding windows (", num.bins, " bins, ", num.steps, " steps)", sep=""), axes=FALSE)
		axis(2)
		for (i in 1:num.steps) {
			for (j in 1:num.bins) {
				if (j == num.bins & i > 1) {
					next
				}
				else {
					rect(xleft=i-1, xright=i, ybottom=bounds[j, i], ytop=bounds[j+1, i], col=colors[num.total-num.steps*(j-1)-i+1], border="white", lwd=2, lty="solid")
				}
			}
		}
		d <- density(x)
		plot(c(min(d$x), max(d$x)), c(0, num.total), type="p", cex=0, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ...)
		axis(1)
		for (step in 1:num.steps) {
			for ( i in 1:num.bins) {
				bin.current <- (i-1)*num.steps+step
				if (bin.current > num.total) {
					next
				}
				if (bounds[i, step] > bounds[i+1, step]) {
					bounds[i+1, step] <- bounds[i, step]
				}
				bin <- which((y >= bounds[i, step]) & (y <= bounds[i+1, step]))
				if (length(bin) < 2) {
					warning("Skipping bin ", bin.current, " (step ", step, ", bin ", bin, "): need at least 2 data points")
					next
				}
				d <- density(x[bin])
				lines(d$x, 0.5*(bin.current)+0.25*(1+(bin.current)/num.total)*num.total*d$y/max(d$y), col=colors[num.total-bin.current+1], lty="solid")
			}
		}
		if (show.avg) {
			d <- density(x)
			lines(d$x, 0.25*num.total+0.5*num.total*d$y/max(d$y), col="gray", lwd=3, lty="solid")
			lines(d$x, 0.25*num.total+0.5*num.total*d$y/max(d$y), col="black", lwd=1, lty="solid")
		}
	}
)

