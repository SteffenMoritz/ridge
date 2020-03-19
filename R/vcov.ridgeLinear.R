## vcov method for ridgeLinear objects

#' @rdname vcov
#' @export
#' @importFrom stats coef model.frame model.matrix .checkMFClasses
vcov.ridgeLinear <- function (object, ...) 
{
  # TODO(dan): check we didn't get any arguments we are unprepared for

  # get back the original data
  data <- model.frame(object)

  # For now, only works if we have an intercept
  stopifnot(object$Inter == 1)
  # drop "(Intercept)" from names
  coef_names <- attr(coef(object), "names")[-1]
  # get data, with intercept term
  X <- model.matrix(object, data=data)
  # make y the original data:
  y_orig <- with(object, y+ym)
  # See also equation 3.44 of "elements of statistical learning", hastie, Tibshirani, Friedman
  # https://web.stanford.edu/~hastie/ElemStatLearn/
  # Inverting this matrix may not be the most numerically stable.
  inv_mat = solve(t(X) %*% X + object$lambda * diag(ncol(X)))
  betaHat <- inv_mat %*% t(X) %*% y_orig
  # sum squared residuals divided by degrees of freedom
  # degrees of freedom: # values minus coefficients and intercept term
  sigma2 <- sum( (y_orig - predict(object, data=data))^2 ) /
    (length(y_orig) - length(coef_names) - 1)
  var_betaHat <- sigma2 * inv_mat
  return(var_betaHat)
}
