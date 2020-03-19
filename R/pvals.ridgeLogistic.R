## computing pvals for ridgeLogistic models

#' @rdname pvals
#' @export
#' @importFrom stats pnorm
pvals.ridgeLogistic <- function(x, ...)
{
  automatic <- x$automatic
  chosen.nPCs <- x$chosen.nPCs
  max.nPCs <- x$max.nPCs
  isScaled <- x$isScaled
  B <- x$coef
  Inter <- x$Inter
  if(Inter)
    {
      X <- cbind(1,x$x)
    } else {
      X <- x$x
    }
  ## lambda may be a vector
  lambda <- x$lambda
  ## B may be a matrix
  xb <- apply(B, 2, function(x){X %*% x})
  expXB <- exp(xb)
  p <- expXB / (1 + expXB)
  W <- vector("list", length = ncol(p))
  for(i in seq(ncol(p)))
    {
      W[[i]] <- diag(p[,i]*(1-p[,i]),length(p[,i]),length(p[,i]))
    }
  KI <- lapply(lambda, function(x){diag(2 * x, dim(X)[2], dim(X)[2])})
  if(Inter)
    {
      for(i in seq(length(lambda)))
      {
        KI[[i]][1,1] <- 0
      }
    }
  computeV <- function(W, KI)
    {
      V <- solve(t(X)%*%W%*%X+KI) %*% (t(X)%*%W%*%X) %*% solve(t(X)%*%W%*%X+KI)
      return(V)
    }
  V <- mapply("computeV", W, KI, SIMPLIFY = FALSE)
  se <- sapply(V, function(x){sqrt(diag(x))})
  tstat <-B/se
  pval <- 2*(1 - pnorm(abs(tstat)))
  if(Inter)
    {
      B <- B[-1, ]
      se <- se[-1, ]
      tstat <- tstat[-1, ]
      pval <- pval[-1, ]
    }
  res <- list(coef = cbind(B), se = cbind(se), tstat = cbind(tstat), pval = cbind(pval), isScaled = isScaled, automatic = automatic, lambda = lambda, chosen.nPCs = chosen.nPCs, max.nPCs = max.nPCs)
  class(res) <- "pvalsRidgeLogistic"
  res
}
