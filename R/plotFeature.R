##DEFINE GENERIC FUNCTION plotFeature()
setGeneric("plotFeature", 
	function(x, y, ...) {
		standardGeneric("plotFeature")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("plotFeature", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("argument 'x' is missing with no default")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotFeature", 
	signature=c("ExpressionSet", "missing"), 
	function(x, y, feature="SIZE", ...) {
		callGeneric(x, getFeatures(x, feature), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("plotFeature", 
	signature=c("ExpressionSet", "vector"), 
	function(x, y, element.x="exprs", element.y="exprs2", sample=1, which.random=NULL, feature.random="TYPE", random.flag="RAND", verbose=FALSE, ...) {
		if (!validObject(x)) {
			stop("not a valid ExpressionSet object")
		}
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Fetching element ", element.x, " for sample ", sample, " from ExpressionSet ...")
		}
		x1 <- getSamples(x, sample, element=element.x)
		if (verbose) {
			cat("FINISHED\nFetching element ", element.y, " for sample ", sample, " from ExpressionSet ...")
		}
		x2 <- getSamples(x, sample, element=element.y)
		if (verbose) {
			cat("FINISHED\n")
		}
		if (is.null(which.random)) {
			if (verbose) {
				cat("Selecting random probes ... ")
			}
			which.random <- getFeatures(x, feature.random)
			if (verbose) {
				cat("FINISHED (", length(which.random), " identified)\nCalling plot() subroutine ... ",sep="")
			}
			callGeneric(cbind(x1, x2), y, which.random=which(which.random == random.flag), verbose=FALSE, ...)
		}
		else {
			if (verbose) {
				cat("Calling plot() subroutine ... ")
			}
			callGeneric(cbind(x1, x2), y, which.random=which.random, verbose=FALSE, ...)
		}
		if (verbose) {
			cat("FINISHED (", (proc.time()["elapsed"]-start), "s elapsed)\n", sep="")
		}
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("plotFeature", 
	signature=c("matrix", "vector"), 
	function(x, y, na.rm=TRUE, limit=10000, which.random=NULL, cutoff=NULL, cutoff2=NULL, main="", xlab="Fragment size (bp)", ylab="log(MspI)", ylab2="log(HpaII)", verbose=FALSE, data.color="black", fail.color="yellow", rand.color="blue", cutoff.color="red", cex=0.2, ...) {
		par.initial <- par()
		if (dim(x)[2] != 2) {
			stop("'x' must contain two columns of data")
		}
		y <- as.vector(y)
		N <- length(y)
		if (dim(x)[1] != N) {
			stop("'x' and 'y' lengths do not match")
		}
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Processing random data ... ")
		}
		which.random <- as.integer(as.vector(which.random))
		which.random <- which.random[which(!is.na(which.random) & (which.random >= 1) & (which.random <= N))]
		random.x <- x[which.random, ]
		random.y <- y[which.random]
		not.random <- which(!(1:N %in% which.random))
		x <- x[not.random, ]
		y <- y[not.random]
		rm(which.random, not.random)
		if (any(is.na(x[, 1]) | is.na(x[, 2]) | is.na(y))) {
			if (na.rm) {
				not.na <- which(!(is.na(x[, 1]) | is.na(x[, 2]) | is.na(y)))
				x <- x[not.na, ]
				y <- y[not.na]
			}
			else {
				stop("NA values where data expected")
			}					
		}
		if (any(is.na(random.x[, 1]) | is.na(random.x[, 2]))) {
			if (na.rm) {
				not.na <- which(!(is.na(random.x[, 1]) | is.na(random.x[, 2])))
				random.x <- random.x[not.na, ]
				random.y <- random.y[not.na]
			}
			else {
				stop("NA values where data expected")
			}					
		}
		N <- length(y)
		if (N < 1) {
			stop("need at least 1 data point to plot")
		}
		if (verbose) {
			cat("FINISHED\nCalculating cutoffs ... ")
		}
		if (is.null(cutoff)) {
			cutoff <- median(random.x[, 1], na.rm=TRUE)+2.5*mad(random.x[, 1], na.rm=TRUE)
		}
		if (is.null(cutoff2)) {
			cutoff2 <- median(random.x[, 2], na.rm=TRUE)+2.5*mad(random.x[, 2], na.rm=TRUE)
		}
		if (verbose) {
			cat("FINISHED\nGenerating plot of features ... ")
		}
		par(mfrow=c(3, 2))
		layout(cbind(c(1, 3, 5), c(2, 4, 6)), widths=rep(c(2, 6), 3), heights=rep(3, 6))
		if (is.null(limit)) {
			limit <- N	
		}
		if (limit > N) {
			limit <- N	
		}
		x.sample <- sample(1:N, size=limit)
		fail.sample <- x.sample[which(x[x.sample, 1] <= cutoff & x[x.sample, 2] <= cutoff2)]
		good <- which((x[, 1] > cutoff) | (x[, 2] > cutoff2))
		fail <- which((x[, 1] <= cutoff) & (x[, 2] <= cutoff2))
		d <- density(x[, 1])
		ylim <- range(d$x)
		plot(max(d$y)-d$y, d$x, type="l", col="white", main="", ylab=ylab, xlab="Density")
		if (length(good) > 2) {
			d1 <- density(x[good, 1])
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=data.color)
		}
		if (length(fail) > 2) {
			d1 <- density(x[fail, 1])
			points(0.5*(2*max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=fail.color)
		}
		if (dim(random.x)[1] > 2) {
			d1 <- density(random.x[, 1], na.rm=TRUE)
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x , type="l", col=rand.color)
		}
		abline(h=cutoff, col=cutoff.color)
		plot(y[x.sample], x[x.sample, 1], type="p", pch=20, cex=cex, main=main, xlab=xlab, ylab="", ylim=ylim, col=data.color, ...)
		points(y[fail.sample], x[fail.sample, 1], type="p", pch=20, cex=cex, col=fail.color)
		abline(h=cutoff, col=cutoff.color)
		d <- density(x[, 2])
		ylim <- range(d$x)
		plot(max(d$y)-d$y, d$x, type="l", col="white", main="", ylab=ylab2, xlab="Density")
		if (length(good) > 2) {
			d1 <- density(x[good, 2])
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=data.color)
		}
		if (length(fail) > 2) {
			d1 <- density(x[fail, 2])
			points(0.5*(2*max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=fail.color)
		}
		if (dim(random.x)[1] > 2) {
			d1 <- density(random.x[, 2], na.rm=TRUE)
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x , type="l", col=rand.color)
		}
		abline(h=cutoff2, col=cutoff.color)
		plot(y[x.sample], x[x.sample, 2], type="p", pch=20, cex=cex, main=main, xlab=xlab, ylab="", ylim=ylim, col=data.color, ...)
		points(y[fail.sample], x[fail.sample, 2], type="p", pch=20, cex=cex, col=fail.color)
		abline(h=cutoff2, col=cutoff.color)
		d <- density(x[, 2] - x[, 1])
		ylim <- range(d$x)
		plot(max(d$y)-d$y, d$x, type="l", col="white", main="", ylab=paste(ylab2, "-", ylab, sep=" "), xlab="Density")
		if (length(good) > 2) {
			d1 <- density(x[good, 2] - x[good, 1])
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=data.color)
		}
		if (length(fail) > 2) {
			d1 <- density(x[fail, 2] - x[fail, 1])
			points(0.5*(2*max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x, type="l", col=fail.color)
		}
		if (dim(random.x)[1] > 2) {
			d1 <- density(random.x[, 2] - random.x[, 1], na.rm=TRUE)
			points((max(d1$y)-d1$y)*(max(d$y)/max(d1$y)), d1$x , type="l", col=rand.color)
		}
		abline(h=mean(random.x[, 2] - random.x[, 1], na.rm=TRUE), col=cutoff.color)
		plot(y[x.sample], x[x.sample, 2] - x[x.sample, 1], type="p", pch=20, cex=cex, main=main, xlab=xlab, ylab="", ylim=ylim, col=data.color, ...)
		points(y[fail.sample], x[fail.sample, 2] - x[fail.sample, 1], type="p", pch=20, cex=cex, col=fail.color)
		abline(h=mean(random.x[, 2] - random.x[, 1], na.rm=TRUE), col=cutoff.color)
		par(par.initial[which(!names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "gamma", "page"))])
		if (verbose) {
			cat("FINISHED (", (proc.time()["elapsed"]-start), "s elapsed)\n", sep="")
		}
	}
)
