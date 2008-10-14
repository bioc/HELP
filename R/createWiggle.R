##DEFINE GENERIC FUNCTION createWiggle()
setGeneric("createWiggle", 
	function(x, y, ...) {
		standardGeneric("createWiggle")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("createWiggle", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		stop("argument 'x' is missing with no default")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("createWiggle", 
	signature=c("ExpressionSet", "missing"), 
	function(x, y, element="exprs", feature.chr="CHR", feature.start="START", feature.stop="STOP", ...) {
		if (!validObject(x)) {
			stop("'x' not a valid ExpressionSet object")
		}
		if (!(element %in% assayDataElementNames(x))) {
			stop("'element' not a valid element of ExpressionSet")
		}
		callGeneric(assayDataElement(x, element), getFeatures(x, c(feature.chr, feature.start, feature.stop)), ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("createWiggle", 
	signature=c("ExpressionSet", "matrix"), 
	function(x, y, element="exprs", ...) {
		if (!validObject(x)) {
			stop("'x' not a valid ExpressionSet object")
		}
		if (!(element %in% assayDataElementNames(x))) {
			stop("'element' not a valid element of ExpressionSet")
		}
		callGeneric(assayDataElement(x, element), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("createWiggle", 
	signature=c("vector", "matrix"), 
	function(x, y, ...) {
		callGeneric(as.matrix(x), y, ...)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "matrix"
setMethod("createWiggle", 
	signature=c("matrix", "matrix"), 
	function(x, y, samples=NULL, na.rm=TRUE, colors=NULL, file=NULL, append=FALSE, sep="\t", ...) {
		x <- getSamples(as.matrix(x), samples, ...)
		N <- dim(x)[2]
		if (N < 1) {
			stop("no data specified")
		}
		if (is.null(colnames(x))) {
			colnames(x) <- 1:N	
		}
		y <- as.matrix(y)
		if (dim(y)[2] != 3) {
			stop("'y' must contain 3 columns of values")
		}
		if (dim(x)[1] != dim(y)[1]) {
			stop("'x' and 'y' lengths differ")
		}
		if (is.null(file)) {
			return()
		}
		if (any(is.na(x))) {
			if (na.rm) {
				not.na <- which(apply(!is.na(x), 1, any))
				x <- x[not.na, ]
				y <- y[not.na, ]
			}
			else {
				stop("'x' contains missing values")
			}
		}
		if (any(is.na(y))) {
			if (na.rm) {
				not.na <- which(apply(!is.na(y), 1, any))
				x <- x[not.na, ]
				y <- y[not.na, ]
			}
			else {
				stop("'y' contains missing values")
			}
		}
		y.order <- order(y[, 1], y[, 2])
		y <- y[y.order, ]
		x <- x[y.order, ]
		scipen <- getOption("scipen")
		options(scipen=10)
		if (is.null(colors)) {
			colors <- rainbow(n=N, start=0, end=2/3)
		}
		if (length(colors) != N) {
			colors <- rainbow(n=N, start=0, end=2/3)
		}
		for (i in 1:length(colors)) {
			colors[i] <- paste(col2rgb(colors[i]), collapse=",")
		}
		if (!append) {
			cat("browser full altGraph\n", file=file)
		}
		for (i in 1:N) {
			cat("track type=wiggle_0 name=\"", colnames(x)[i], "\" description=\"HELP\" visibility=full autoScale=off color=", colors[i], " altColor=", colors[i], " priority=", i, " yLineOnOff=on yLineMark=0\n", file=file, append=TRUE, sep="")
			apply(cbind(y[, 1:3], x[, i]), 1, cat, "\n", file=file, append=TRUE, sep=sep)
		}
		options(scipen=scipen)
	}
)
