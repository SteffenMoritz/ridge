## plot the ridge trace
## for ridgeLogistic objects

#' @rdname plot
#' @export
#' @importFrom stats coef
#' @importFrom graphics plot lines abline
#' @importFrom grDevices rainbow
plot.ridgeLogistic <- function(x, y = NULL, ...)
  {
    Inter <- x$Inter
    lambda <- x$lambda
    automatic <- x$automatic
    nPCs <- x$max.nPCs
    coefs <- rbind(coef(x, all.coef = TRUE))
    if(Inter)
      {
        coefs <- coefs[,-Inter]
      }
    ## x is a ridgeLinear object
    if(length(lambda) == 1)
      {
        col.vector <- rainbow(length(coefs))
        plot(x = rep(lambda, length(coefs)), y = coefs, xlab = "lambda", ylab = "coefficient", col = col.vector, pch = 19)
      } else {
        col.vector <- rainbow(ncol(coefs))
        if (automatic) {
          chosen.nPCs <- x$chosen.nPCs
          plot(x = seq(nPCs), y = coefs[,1], ylim = range(coefs), xlab = "nPCs", ylab = "coefficient", col = col.vector[1], type = "l", main = "ridge trace")
          for(i in 2:ncol(coefs))
            {
              lines(x = seq(nPCs), y = coefs[,i], col = col.vector[i])
            }
          abline(v = chosen.nPCs, lty = 2)
        } else {
          plot(x = lambda, y = coefs[,1], xlim=range(lambda), ylim = range(coefs), xlab = "lambda", ylab = "coefficient", col = col.vector[1], type = "l", main = "ridge trace")
          for(i in 2:ncol(coefs))
            {
              lines(x = lambda, y = coefs[,i], col = col.vector[i])
            }
        }
      }
  }
