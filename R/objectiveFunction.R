objectiveFunction <- function(B, X, y, k, intercept = TRUE)
  {
    #Numeric bits for objective function
    #compute p
    XB <- X%*%B
    expXB <- exp(XB)
    p <- expXB/(1+expXB)
    p[is.infinite(expXB)] <- 1
    top <- ifelse(y,log(p),log(1-p))
    if(intercept)
      {
        objBeta <- sum(top) - k*crossprod(B[-1],B[-1])
      } else {
        objBeta <- sum(top) - k*crossprod(B,B)
      }
    return(objBeta)
  }
