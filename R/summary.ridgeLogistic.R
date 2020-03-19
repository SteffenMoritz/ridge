## summary function for ridgeLogistic object

#' @rdname summary
#' @export
#' @importFrom stats pnorm
summary.ridgeLogistic <- function(object, all.coef = FALSE, ...)
  {
    res <- vector("list")
    isScaled <- object$isScaled
    Inter <- object$Inter
    res$automatic <- object$automatic
    res$call <- object$call
    res$lambda <- object$lambda
    pvalues <- pvals(object)
    summaries <- vector("list", length(res$lambda))
    res$all.coef = all.coef
    coefs <- rbind(coef(object, all.coef = TRUE))
    if(res$automatic)
      {
        res$chosen.nPCs <- object$chosen.nPCs
      }
    for(i in seq(length(res$lambda)))
      {
        summary <- vector("list")
        if(Inter)
          {
            if(isScaled) {
              ## Both intercept and scaled
              summary$coefficients <- cbind(coefs[i,], c(NA, pvalues$coef[,i]), c(NA, pvalues$se[,i]), c(NA, pvalues$tstat[,i]), c(NA, pvalues$pval[,i]))
              dimnames(summary$coefficients) <- list(c("(Intercept)", names(object$xm)), c("Estimate", "Scaled estimate", "Std. Error (scaled)", "t value (scaled)", "Pr(>|t|)"))
            } else {
              ## Intercept, not scaled
              summary$coefficients <- cbind(coefs[i,], c(NA, pvalues$se[,i]), c(NA, pvalues$tstat[,i]), c(NA, pvalues$pval[,i]))
              dimnames(summary$coefficients) <- list(c("(Intercept)", colnames(object$x)), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
            }
          } else {
            if(isScaled) {
              summary$coefficients <- cbind(coefs[i,], pvalues$coef[,i], pvalues$se[,i], pvalues$tstat[,i], pvalues$pval[,i])
              dimnames(summary$coefficients) <- list(colnames(object$x), c("Estimate", "Scaled estimate", 
                                                                       "Std. Error (scaled)", "t value (scaled)", "Pr(>|t|)"))
            } else {
              summary$coefficients <- cbind(coefs[i,], pvalues$se[,i], pvalues$tstat[,i], pvalues$pval[,i])
              dimnames(summary$coefficients) <- list(colnames(object$x), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
            }
          }
        ## ## df: degrees of freedom: model, variance, residual (may as well return all three)
        summary$df <- object$df[i,]
        ## nPCs: number of principal components used to choose lambda
        if(!is.null(object$max.nPCs))
          {
            summary$nPCs <- i
          }
        summary$lambda <- object$lambda[i]
        summaries[[i]] <- summary
        names(summaries)[[i]] <- paste("summary", i, sep = "")
        rm(summary)
      } ## ends for(i in seq(length(res$lambfa)))
    res$summaries <- summaries
    ## Make an object of class summary.ridgeLogistic
    class(res) <- "summary.ridgeLogistic"
    ## Call its print method (print.summary.ridgeLogistic)
    res
  }
