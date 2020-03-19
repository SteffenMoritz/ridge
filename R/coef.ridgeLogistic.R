## coef method for ridgeLogistic objects

#' @rdname coef
#' @export
coef.ridgeLogistic <- function (object, all.coef = FALSE, ...) 
{
  if (object$Inter) {
    Intercept <- apply(object$coef, 2, function(x) {x[1] - sum(object$xm * x[-1] / object$scales)})
    scaledcoefs <- rbind(Intercept, apply(object$coef, 2, function(x){x[-1] / object$scales}))
    rownames(scaledcoefs)[1] <- "(Intercept)"
   } else {
    scaledcoefs <- apply(object$coef, 2, function(x){x/object$scales})
  }
  if(object$automatic && all.coef == FALSE)
    {
      scaledcoefs <- scaledcoefs[,object$chosen.nPCs]
    }
  scaledcoefs <- t(scaledcoefs)
  drop(scaledcoefs)
}
