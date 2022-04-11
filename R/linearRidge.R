## R function to fit the linear ridge regression model

#' @export
#' @importFrom stats .getXlevels predict model.response model.matrix
linearRidge <- function(formula, data, lambda = "automatic",
                        nPCs = NULL, scaling = c("corrForm", "scale", "none"), ...)
  {
    ## Check arguments
    if(lambda != "automatic" && !is.null(nPCs))
      {
        stop(gettextf("you cannot specify both lambda and nPCs\n"))
      } else if(lambda == "automatic" && !is.null(nPCs)) {
        lambda <- NULL
      }
    automatic <- FALSE
    cl <- match.call()
    m <- match.call(expand.dots = FALSE)
    scaling <- match.arg(scaling)
    if((lambda == "automatic" && !is.null(lambda)) && scaling != "corrForm")
      {
        warning("lambda is chosen automatically so scaling is set to \"corrForm\" ")
        scaling <- "corrForm"
        cl$scaling <- scaling
      }
    if(is.null(lambda) && scaling != "corrForm")
      {
        warning("lambda is chosen based on number of components so scaling is set to \"corrForm\"")
        scaling <- "corrForm"
        cl$scaling <- scaling
      }
    if(scaling == "none")
      {
        isScaled <- FALSE
        corrForm <- FALSE
        standardize <- FALSE
      } else if (scaling == "corrForm") {
        isScaled <- TRUE
        corrForm <- TRUE
        standardize <- FALSE
      } else if (scaling == "scale") {
        isScaled <- TRUE
        corrForm <- FALSE
        standardize <- TRUE
      }
    m$model <- m$allLambdas <- m$nPCs <- m$... <- m$lambda <- m$scaling <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")

    ## Extract the response
    Y <- model.response(m)
    ## Construct the design matrix
    X <- model.matrix(Terms, m)
    contrasts <- attr(X, "contrasts")
    xlevels <- .getXlevels(Terms, m)
    ## get the dimensions of X in terms of n and p
    n <- nrow(X)
    p <- ncol(X)
    if(!is.null(nPCs) && nPCs > min(n, p))
      stop(gettextf("You specified %d PCs which is greater than the maximum number of PCs in the data\n", nPCs))
    ## Handle the intercept
    if (Inter <- attr(Terms, "intercept"))
      {
        Xm <- colMeans(X[, -Inter])
        Ym <- mean(Y)
        p <- p - 1
        ## Subtract the means from X
        X <- X[, -Inter] - rep(Xm, rep(n, p))
        ## Subtract the mean from Y
        Y <- Y - Ym
      } else {
        Xm <- colMeans(X)
        Ym <- mean(Y)
        ## Subtract the means from X
        X <- X - rep(Xm, rep(n, p))
        ## Subtract the mean from y
        Y <- Y - Ym
      }
	## Calculate the scales
     if(corrForm)
       {
         Xscale <- drop(rep(1/(n - 1), n) %*% apply(X, 2, function(x){x - mean(x)})^2)^0.5 * sqrt(nrow(X) - 1)
       } else if (standardize) {
         Xscale <- drop(rep(1/(n - 1), n) %*% apply(X, 2, function(x){x - mean(x)})^2)^0.5
       } else {
         Xscale <- drop(rep(1, p))
         names(Xscale) <- colnames(X)
       }
     X <- X/rep(Xscale, rep(n, p))
     Xs <- svd(X)
     Q <- Xs$v
     ## Make the principal components
     Z <- X %*% Q
     Lambda <- Xs$d^2
     if(!is.null(lambda) && lambda == "automatic")
       {
         automatic <- TRUE
         if(is.null(nPCs))
           {
             propVar <- cumsum(Lambda) / sum(Lambda) * 100
             # if there are values for propVar below 90
             if (length(which.max(propVar[propVar < 90])) > 0)
             {
               ifelse((length(propVar[propVar >= 90]) > 0), (nPCs <- which.max(propVar[propVar < 90]) + 1), (nPCs <- which.max(propVar[propVar < 90])))
             }
             else
             {
               nPCs <- 1
             }
             
           }
       }
    ## Compute ahat
    ahat <- diag(1 / Lambda) %*% t(Z) %*% Y
    if(!is.null(lambda) && lambda == "automatic" && !is.null(nPCs))
      {
        ks.vector <- sig2hat.vector <- vec.df <- numeric(nPCs)
        flag <- TRUE
        P <- 0
        while((P < nPCs) && flag)
          {
            P <- P + 1
            ## compute sig2hatP
            sig2hat <- ifelse(P == 1,
                              as.numeric(crossprod(Y - (Z[,1]) * ahat[1]) / (n - 1)),
                              as.numeric(crossprod(Y - (Z[,1:P]) %*% ahat[1:P]) / (n - P))
                              )
            ## compute ahatsum
            ahatsum <- ifelse(P == 1,
                              ahat[1]^2,
                              sum(ahat[1:P]^2)
                              )
            ## compute kHKB
            ks.vector[P] <- P * sig2hat / ahatsum
            if(is.finite(ks.vector[P]))
              {
                vec.df[ P ] <- sum(Lambda^2 / (Lambda + ks.vector[P])^2)
              }
            if(!is.finite(ks.vector[P]))
              {
                flag <- FALSE
                ## make everything the correct dimensions
                ks.vector <- ks.vector[1:(P - 1)]
              }
          } ## Ends while loop
        ## Choose best lambda
        nPCs.dof <- which.min(abs(vec.df - seq(nPCs)))
        ## Vector of lambdas
        lambda <- ks.vector
        ## The number of components
        chosen.nPCs <- nPCs.dof
        max.nPCs <- nPCs
      } else if (!is.null(nPCs))
        {
          P <- nPCs
          sig2hat <- ifelse(P == 1,
                            as.numeric(crossprod(Y - (Z[,1]) * ahat[1]) / (n - 1)),
                            as.numeric(crossprod(Y - (Z[,1:P]) %*% ahat[1:P]) / (n - P))
                            )
          ahatsum <- ifelse(P == 1,
                            ahat[1]^2,
                            sum(ahat[1:P]^2)
                            )
          ## compute lambda
          lambda <- P * sig2hat / ahatsum
          chosen.nPCs <- nPCs
        }
    ## compute coef as a matrix
    aridge <- lapply(lambda, function(x) {ahat * Lambda / (Lambda + x)})
    coef <- lapply(aridge, function(x) {Q %*% x})
    ## compute df as a matrix
    df <- lapply(lambda, function(x) {c(sum(Lambda / (Lambda + x)), sum(Lambda^2 / (Lambda + x)^2), sum(Lambda * (Lambda + 2* x) / (Lambda + x)^2))})
    coef <- do.call(cbind, coef)
    rownames(coef) <- colnames(X)
    ##
    if(!is.null(nPCs))
      {
        if(length(lambda) == 1)
          {
            colnames(coef) <- paste("nPCs", chosen.nPCs, sep = "")
          } else {
            colnames(coef) <- paste("nPCs", seq(max.nPCs), sep = "")
          }
      } else {
        colnames(coef) <- paste("lambda=", lambda, sep = "")
      }
    ##
    df <- do.call(rbind, df)
    ## This line needs fixing
    if(!is.null(nPCs))
      {
        if(length(lambda) == 1)
          {
            rownames(df) <- paste("nPCs", chosen.nPCs, sep = "")
          } else {
            rownames(df) <- paste("nPCs", seq(max.nPCs), sep = "")
          }
      } else {
        rownames(df) <- paste("lambda=", lambda, sep = "")
      }
    ##
    colnames(df) <- c("model", "variance", "residual")
    res <- list(automatic = automatic, call = cl, coef = cbind(drop(coef)), df = df,
                Inter = Inter, isScaled = isScaled, lambda = lambda, scales = Xscale,
                terms = Terms, x = X, xm = Xm, y = Y, ym = Ym, model_frame = m,
                contrasts = contrasts, xlevels = xlevels)
    ## This line needs fixing
    if(!is.null(nPCs))
      {
        if(automatic)
          {
            res$max.nPCs <- max.nPCs
          }
        res$chosen.nPCs <- chosen.nPCs
      }
    class(res) <- "ridgeLinear"
    res
  }
