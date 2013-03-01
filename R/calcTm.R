##DEFINE GENERIC FUNCTION calcTm()
setGeneric("calcTm", 
	function(x, ...) {
		data(base.stacking.thermodynamics)
		standardGeneric("calcTm")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("calcTm", 
	signature=c("missing"), 
	function(x, ...) {
		return(NULL)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "NULL"
setMethod("calcTm", 
	signature=c("NULL"), 
	function(x, ...) {
		return(NULL)
	}
)


##DEFINE METHOD TO HANDLE CLASS: "character"
setMethod("calcTm", 
	signature=c("character"), 
	function(x, strand1.concentration=2e-07, strand2.concentration=2e-07, method="nearest-neighbor", ...) {
		x <- toupper(x)
		seq.lengths <- nchar(x)
		sequence <- rep(NA, length(x))
		chars.all <- strsplit(x, split="")
		for (curr in 1:length(chars.all)) {
			chars <- chars.all[[curr]]
			if (any(!(chars %in% c("A", "T", "C", "G")))) {
				warning("Unable to calculate Tm: unrecognized characters")
				next	
			}
			seq.length <- seq.lengths[curr]
			if (seq.length <= 2) {
				warning("Unable to calculate Tm: sequence too short")
				next			
			}
			dH <- 0
			dS <- -10.8
			if ((chars[1] == "A") | (chars[1] == "T")) {
				dS <- dS+4.1
				dH <- dH+2.3
			}
			else {
				dS <- dS-2.8
				dH <- dH+0.1
			}
			if ((chars[seq.length] == "A") | (chars[seq.length] == "T")) {
				dS <- dS+4.1
				dH <- dH+2.3
			}
			else {
				dS <- dS-2.8
				dH <- dH+0.1
			}
			if (all(chars == rev(chars))) {
				dS <- dS-1.4
			}
			for (i in 1:(seq.length-1)) {
				seq.dinucleotide <- substr(x[curr],i,i+1)
				dH <- dH+base.stacking.thermodynamics[seq.dinucleotide,1]
				dS <- dS+base.stacking.thermodynamics[seq.dinucleotide,2]
			}
			sequence[curr] <- as.double(1000*dH/(dS+1.987*log(max(strand1.concentration, strand2.concentration)-min(strand1.concentration, strand2.concentration)/2)) - 273.15)
		}
		return(as.numeric(sequence))
	}
)



##DEFINE METHOD TO HANDLE CLASS: "character"
setMethod("calcTm", 
	signature=c("ExpressionSet"), 
	function(x, ...) {
		if (!validObject(x)) {
			stop("'x' not a valid ExpressionSet object")
		}
		if (!"SEQUENCE" %in% varLabels(featureData(x))) {
			stop("featureData must contain 'SEQUENCE' information")
		}
		if ("TM" %in% varLabels(featureData(x))) {
			stop("featureData already contains Tm information")
		}
		sequence <- as.character(pData(featureData(x))[, "SEQUENCE"])
		tm <- callGeneric(sequence, ...)
		tm <- data.frame(TM=as.numeric(tm))
		featureData(x) <- combine(featureData(x), new("AnnotatedDataFrame", data=tm, dimLabels=c("featureNames", "featureColumns")))
		return(x)
	}
)
