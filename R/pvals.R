## Generic method for computing pvalues on ridgeLinear or ridgeLogistic models

pvals <- function(x, ...)
	UseMethod("pvals")
	
