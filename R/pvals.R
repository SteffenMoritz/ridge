## Generic method for computing pvalues on ridgeLinear or ridgeLogistic models

#' @export
pvals <- function(x, ...)
	UseMethod("pvals")
	
