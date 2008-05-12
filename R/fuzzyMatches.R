##DEFINE GENERIC FUNCTION fuzzyMatches()
setGeneric("fuzzyMatches", 
	function(x, y, ...) {
		standardGeneric("fuzzyMatches")
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("fuzzyMatches", 
	signature=c("missing", "missing"), 
	function(x, y, ...) {
		return()
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("fuzzyMatches", 
	signature=c("vector", "missing"), 
	function(x, y, ...) {
		return()
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("fuzzyMatches", 
	signature=c("vector", "NULL"), 
	function(x, y, ...) {
		return()
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("fuzzyMatches", 
	signature=c("NULL", "vector"), 
	function(x, y, ...) {
		return()
	}
)

##DEFINE METHOD TO HANDLE CLASS: "vector"
setMethod("fuzzyMatches", 
	signature=c("vector", "vector"), 
	function(x, y, strict=TRUE, keep=TRUE, alias=TRUE, match.all="*", nomatch=0, na.rm=TRUE, verbose=FALSE, ...) {
		if (verbose) {
			start <- proc.time()["elapsed"]
			cat("Matching (", paste(x, collapse=", "), ") against ", paste(y, collapse=", "), ") ... ", sep="")
		}
		x <- unique(x)
		y <- unique(y)
		if (na.rm) {
			x <- x[which(!is.na(x))]	
			y <- y[which(!is.na(y))]	
		}
		if (!is.null(match.all) & any(x == match.all) & !strict) {
			x.index <- which(x == match.all)
			if (x.index == 1) {
				if (length(x) > 1) {
					x <- c(y, x[2:length(x)])	
				}
				else {
					x <- y	
				}
			}
			else {
				if (x.index == length(x)) {
					x <- c(x[1:(x.index - 1)], y)
				}
				else {
					x <- c(x[1:(x.index - 1)], y, x[(x.index+1):length(x)])	
				}			
			}
		}
		matches <- match(x, y)
		if (verbose) {
			cat(length(which(!is.na(matches))), " exact, ")
		}
		x.index <- x
		x.index[which(!is.na(matches))] <- matches[which(!is.na(matches))]
		if (strict) {
			x <- x[which(!is.na(matches))]	
			x.index <- x.index[which(!is.na(matches))]
		}
		x.index <- as.integer(x.index)
		if (any(is.na(x.index))) {
			warning("unmatched values in 'x': ", paste(x[which(is.na(x.index))], collapse=", "))
		}
		if (verbose) {
			cat(length(which(!is.na(x.index))), " by index, ")
		}
		if (!keep) {
			x <- x[which(!is.na(x.index) & (1 <= x.index & x.index <= length(y)))]	
			x.index <- x.index[which(!is.na(x.index) & (1 <= x.index & x.index <= length(y)))]
		}
		x.index <- x.index[which(is.na(x.index) | (1 <= x.index & x.index <= length(y)))]
		x[which(!is.na(x.index))] <- y[x.index[which(!is.na(x.index))]]
		x <- unique(x)
		if (!alias) {
			matches <- match(x, y)	
			x[which(!is.na(matches))] <- matches[which(!is.na(matches))]
			x[which(is.na(matches))] <- nomatch
		}
		if (length(x) < 1) {
			x <- NULL	
		}
		return(x)
	}
)
