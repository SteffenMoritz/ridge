## coef method for ridgeLinear objects

#' @rdname coef
#' @export
coef.ridgeLinear <- function (object, all.coef = FALSE, ...) 
{
  scaledcoef <- t(as.matrix(object$coef/object$scales))
  if (object$Inter) {
    inter <- object$ym - scaledcoef %*% object$xm
    scaledcoef <- cbind(Intercept = inter, scaledcoef)
    colnames(scaledcoef)[1] <- "(Intercept)"
  }
  if(object$automatic && all.coef == FALSE)
    {
      scaledcoef <- scaledcoef[object$chosen.nPCs,]
    }
  drop(scaledcoef)
}
