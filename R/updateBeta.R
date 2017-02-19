updateBeta <- function(B, X, y, k, intercept = TRUE, doff = FALSE)
  {
    XB <- X%*%B
    expXB <- exp(XB)
    p <- expXB/(1+expXB)
    p[is.infinite(expXB)] <- 1
    W <- diag(as.numeric((p*(1-p))),length(y),length(y))
    WZ <- p*(1-p)*XB+(y-p)    
    kI <- diag(2*k,dim(X)[2],dim(X)[2])
    if(intercept)
      kI[1,1] <- 0 ##Intercept
    updatedBeta <- (solve(t(X)%*%W%*%X+kI))%*%t(X)%*%WZ
    if(doff)
      {
        res = list(updatedBeta = updatedBeta, kI = kI, W = W)
      } else {
        res = updatedBeta
      }
    return(res)
  }
