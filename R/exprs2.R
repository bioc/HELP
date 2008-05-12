##DEFINE GENERIC FUNCTION exprs2()
setGeneric("exprs2", 
	function(object) {
		standardGeneric("exprs2")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("exprs2", 
	signature=c("missing"), 
	function(object) {
		return(NULL)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("exprs2", 
	signature=c("ExpressionSet"), 
	function(object) {
		assayDataElement(object, "exprs2")
	}
)


##DEFINE GENERIC FUNCTION exprs2()
setGeneric("exprs2<-", 
	function(object, value) {
		standardGeneric("exprs2<-")
	}
)


##DEFINE METHOD TO HANDLE CLASS: "missing"
setMethod("exprs2<-", 
	signature=c("ExpressionSet", "missing"), 
	function(object, value) {
		return(object)
	}
)

##DEFINE METHOD TO HANDLE CLASS: "ExpressionSet"
setMethod("exprs2<-", 
	signature=c("ExpressionSet", "matrix"), 
	function(object, value) {
		assayDataElement(object, "exprs2") <- value
		return(object)
	}
)



