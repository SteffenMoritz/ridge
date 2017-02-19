computeRidgeLogistic <- function(X, y, k, intercept = TRUE, doff = FALSE)
  {
    if(is.null(ncol(X)))
      {
        X <- cbind(X)
      }
    if(is.null(colnames(X)))
      {
        colnames(X) <- paste("pred", seq(ncol(X)), sep = "")
      }
    if(intercept)
      {
        X <- cbind(1, X)
        colnames(X)[1] <- "(Intercept)"
      }
    ##Initialize the estimate of B
    B <- numeric(dim(X)[2])
    ##Initial objective function
    objOld <- objectiveFunction(B, X, y, k, intercept)
    ##Updated estimates of beta
    newB <- updateBeta(B, X, y, k, intercept, doff)
    if(doff)
      {
        objNew <- objectiveFunction(newB$updatedBeta, X, y, k, intercept)
      } else {
        objNew <- objectiveFunction(newB, X, y, k, intercept)
      }
    index <- 0
    while(abs(diff(c(objOld,objNew)))>10^-6)
      {
        index <- index+1
        objOld <- objNew
        if(doff)
          {
            B <- newB$updatedBeta
          } else {
            B <- newB
          }
        newB <- updateBeta(B, X, y, k, intercept, doff)
        if(doff){
         objNew <- objectiveFunction(newB$updatedBeta, X, y, k, intercept)
        } else {
          objNew <- objectiveFunction(newB, X, y, k, intercept)
        }
      }
    if(doff)
      {
        B <- newB$updatedBeta
        W <- newB$W
        kI <- newB$kI
      } else {
        B <- newB
      }
    if(doff)
      {
        H <- W %*% X %*% solve(t(X) %*% W %*% X + kI) %*% t(X)
        doff <- c( sum(diag(H)) , sum(diag(H %*% t(H))))
        return(list(B = B, doff = doff))
      }
    else
      return(B)
  }
