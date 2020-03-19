#' @rdname nobs
#' @export
nobs.ridgeLinear <- function(object, ...) {
  length(object$y)
}
