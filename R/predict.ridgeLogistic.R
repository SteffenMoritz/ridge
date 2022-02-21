## predict method for logistic ridge regression models

## fitted.values are the predicted probabilities
## i.e. exp(XB) / (1 + exp(XB))
## linear.predictors are the scores (X %*% B)
## on the scale of the original data

## predict.glm code

#' @rdname predict
#' @export
#' @importFrom stats na.pass terms model.frame delete.response .checkMFClasses model.matrix coef
predict.ridgeLogistic <- function (object, newdata = NULL, type = c("link", "response"), 
    na.action = na.pass, all.coef = FALSE, ...) 
{
  tt <- terms(object)
  type <- match.arg(type) ## Match the type argument
  na.act <- object$na.action ## Get the na.action statement
  object$na.action <- NULL ## Set object$na.action to NULL
  
  if (missing(newdata) || is.null(newdata)) { ## If there is no newdata
    newdata <- object$model_frame
  } 
  
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    mm <- X <- model.matrix(Terms, m)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)

  hasintercept <- attr(tt, "intercept")
  
  ll <- attr(tt, "term.labels")
  
  if(hasintercept)
    mm <- cbind(1, X[,ll])
  else
    mm <- X[,ll]
  B <- coef(object, all.coef = all.coef)
  if(all.coef)
    {
      XB <- apply(B, 1, function(x) {as.matrix(X) %*% x})
    } else {
      XB <- as.matrix(X) %*% B
    }
  expXB <- exp(XB)
  p <- expXB / (1 + expXB)
  pred <- switch(type, link = XB, response = p)
  
  pred
}
