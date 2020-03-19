## computing pvals for ridgeLinear models

#' @rdname pvals
#' @export
#' @importFrom stats pnorm
pvals.ridgeLinear <- function(x, ...)
{
  automatic <- x$automatic
  chosen.nPCs <- x$chosen.nPCs
  max.nPCs <- x$max.nPCs
  isScaled <- x$isScaled
  beta <- x$coef
  names(beta) <- colnames(x$x)
  y <- x$y
  n <- length(y)
  svdX <- svd(x$x)
  U <- svdX$u
  D <- svdX$d
  D2 <- svdX$d^2
  V <- svdX$v
  lambda <- x$lambda
  div <- lapply(lambda, function(x) {D2 + x})
  sig2hat <- apply(rbind(lambda, do.call(cbind, div)), 2, function(x) {as.numeric(crossprod(y - U %*% diag((D2)/(x[-1])) %*% t(U) %*% y)) / (n - sum(D2 * (D2 + 2 * x[1]) / (x[-1]^2)))})
  varmat <- lapply(div, function(x){V %*% diag(D2 / (x^2)) %*% t(V)})
  varmat <- mapply(function(x, y){x * y}, sig2hat, varmat, SIMPLIFY = FALSE)
  se <- lapply(varmat, function(x){sqrt(diag(x))})
  se <- do.call(cbind, se)
  rownames(se) <- rownames(beta)
  colnames(se) <- colnames(beta)
  tstat <- abs(beta / se)
  pval <- 2 * (1 - pnorm(tstat))
  res <- list(coef = beta, se = se, tstat = tstat, pval = pval, isScaled = isScaled, automatic = automatic, lambda = lambda, chosen.nPCs = chosen.nPCs, max.nPCs = max.nPCs)
  class(res) <- "pvalsRidgeLinear"
  res
}
