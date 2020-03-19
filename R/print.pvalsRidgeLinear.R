## print method for pvalsRidgeLinear objects

#' @rdname print
#' @export
#' @importFrom stats printCoefmat
print.pvalsRidgeLinear <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), all.coef = FALSE, ...)
  {
    if(x$automatic && all.coef == FALSE)
      {
        i <- x$chosen.nPCs
        coefs <- cbind(x$coef[,i], x$se[,i], x$tstat[,i], x$pval[,i])
        rownames(coefs) <- rownames(x$coef)
        if(x$isScaled)
          {
            colnames(coefs) <- c("Estimate (scaled)", "Std. Error (scaled)", "t value (scaled)", "Pr(>|t|)")
          } else {
            colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
          }
        cat(gettextf("\n\tlambda %f chosen automatically using %d PCs\n\n", x$lambda[i], x$chosen.nPCs))
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                     na.print = "NA", ...)
      } else {
        for(i in seq(length(x$lambda)))
          {
            coefs <- cbind(x$coef[,i], x$se[,i], x$tstat[,i], x$pval[,i])
            rownames(coefs) <- rownames(x$coef)
            if(x$isScaled)
              {
                colnames(coefs) <- c("Estimate (scaled)", "Std. Error (scaled)", "t value (scaled)", "Pr(>|t|)")
              } else {
                colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
              }
            cat(gettextf("\nlambda %f", x$lambda[i]))
            if(x$automatic && (x$chosen.nPCs == i))
              cat(", chosen automatically")
            if(!is.null(x$max.nPCs))
              cat(gettextf(", computed using %d PCs\n", i))
            cat("\n")
            printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                         na.print = "NA", ...)
          }
      }
    invisible(x)
  }
