#' @rdname nobs
#' @export
#' @importFrom stats nobs
 
nobs.ridgeLinear <- function(object, ...) {
  length(object$y)
}
