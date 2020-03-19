## Plot the pval trace
## For pvalsRidgeLinear objects

#' @rdname plot
#' @export
plot.pvalsRidgeLinear <- function(x, y = NULL, ...)
  {
    lambda <- x$lambda
    automatic <- x$automatic
    nPCs <- x$max.nPCs
    pval <- -1 * log10(x$pval)
    pval[is.infinite(pval)] <- NA
    ## x is a pvalsRidgeLinear object
    if(length(lambda) == 1)
      {
        col.vector <- grDevices::rainbow(length(pval))
        graphics::plot(x = rep(lambda, length(pval)), y = pval, xlab = "lambda", ylab = "-log(10) pvalue", col = col.vector, pch = 19, ...)
      } else {
        col.vector <- grDevices::rainbow(nrow(pval))
        if (automatic) {
          chosen.nPCs <- x$chosen.nPCs
          graphics::plot(x = seq(nPCs), y = pval[1,], ylim = c(0, max(pval, na.rm = TRUE)), xlab = "nPCs", ylab = "-log(10) pvalue", col = col.vector[1], type = "l", main = 
          "pvalue trace", ...)
          for(i in 2:nrow(pval))
            {
              graphics::lines(x = seq(nPCs), y = pval[i,], col = col.vector[i], ...)
            }
          graphics::abline(v = chosen.nPCs, lty = 2)
        } else {
          graphics::plot(x = lambda, y = pval[1,], xlim=range(lambda), ylim = c(0, max(pval, na.rm = TRUE)), xlab = "lambda", ylab = "-log(10) pvalue", col = col.vector[1], 
          type = "l", main = "pvalue trace", ...)
          for(i in 2:nrow(pval))
            {
              graphics::lines(x = lambda, y = pval[i,], col = col.vector[i], ...)
            }
        }
      }
  }
