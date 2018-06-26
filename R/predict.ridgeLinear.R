## predict method for linear ridge regression models

predict.ridgeLinear <- function(object, newdata,  
    na.action = na.pass, all.coef = FALSE, ...)
  {
    tt <- terms(object)
    if (!inherits(object, "ridgeLinear")) 
      warning("calling predict.ridgeLinear(<fake-ridgeLinear-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
      mm <- X <- model.frame(object)
      mmDone <- TRUE
      offset <- object$offset
    }
    else {
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
      mmDone <- FALSE
    }
    hasintercept <- attr(tt, "intercept")
    ll <- attr(tt, "term.labels")
    if(hasintercept)
      mm <- cbind(1, X[,ll,drop=FALSE])
    else
      mm <- X[,ll,drop=FALSE]
    beta <- coef(object, all.coef = all.coef)
    if(all.coef)
      res <- apply(beta, 1, function(x){drop(as.matrix(mm) %*% x)})
    else
      res <- drop(as.matrix(mm) %*% beta)
    res
  }
