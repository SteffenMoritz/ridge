## R function to fit the logistic ridge regression model

#' @export
#' @importFrom stats model.response model.matrix
logisticRidge <- function(formula, data, lambda = "automatic",
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
    m$model <- m$nPCs <- m$... <- m$lambda <- m$scaling <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    
    ## Extract the response
    Y <- model.response(m)
    ## Construct the design matrix
    X <- model.matrix(Terms, m)
    ## get the dimensions of X in terms of n and p
    n <- nrow(X)
    p <- ncol(X)
    if(!is.null(nPCs) && nPCs > min(n, p))
      stop(gettextf("You specified %d PCs which is greater than the maximum number of PCs in the data\n", nPCs))
    ## Handle the intercept
    if (Inter <- attr(Terms, "intercept"))
      {
        Xm <- colMeans(X[, -Inter])
        p <- p - 1
        ## Subtract the means from X
        X <- X[, -Inter] - rep(Xm, rep(n, p))
      }
     else Ym <- Xm <- NA ## Else Ym and Xm are not needed
    ## Because an intercept does not have to be calculated
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
               ifelse((length(propVar[propVar >= 90]) > 0), (max.nPCs <- which.max(propVar[propVar < 90]) + 1), (max.nPCs <- which.max(propVar[propVar < 90])))
             }
             else 
             {
               max.nPCs <- 1
             }
           }
       }
    if(lambda == "automatic" && is.null(nPCs))
      {
        ks.vector <- vec.df <- numeric(max.nPCs)
        flag <- TRUE
        P <- 0
        while(P < max.nPCs && flag)
          {
            P <- P + 1
            tryCatch({
              ks.vector[P] <- P / sum(computeRidgeLogistic(Z[,1:P], Y, 0, intercept = FALSE)^2)
              fittedB <- computeRidgeLogistic(X, Y, k = ks.vector[P], intercept = Inter, doff = TRUE)
              vec.df[P] <- fittedB$doff[2]
            }, error = function(e) {
              flag <- FALSE
            }
                     )
          } ## Ends while
        options("warn" = -1)
        nPCs <- ifelse(!is.infinite(min(which(ks.vector == 0))),
                      min(which(ks.vector == 0)) - 1,
                      max.nPCs)
        options("warn" = 0)
        ks.vector <- ks.vector[1:nPCs]
        vec.df <- vec.df[1:nPCs]
        lambda <- ks.vector
        ## Choose best lambda
        chosen.nPCs <- which.min(abs(vec.df - seq(nPCs) - Inter)) ## The -1 is for the intercept which was included in doff
        max.nPCs <- nPCs
        ## Ends if Lambda == automatic
      } else if (!is.null(nPCs)) {
        P <- nPCs
        tryCatch({
          lambda <- P / sum(computeRidgeLogistic(Z[,1:P], Y, 0, intercept = FALSE)^2)
}, error = function(e) {
              stop(gettextf("Unable to fit logistic ridge model using %d components\n", nPCs))
            })
        chosen.nPCs <- nPCs
      }
    ## Need to make a matrix of fitted B
    ## Don't scale X because X are already scaled and because B on scaled data is needed to compute p-values
    fittedB <- lapply(lambda, function(x){computeRidgeLogistic(X, Y, k = x, intercept = Inter, doff = TRUE)})
    coef <- lapply(fittedB, function(x){x$B})
    df <- lapply(fittedB, function(x){x$doff})
    coef <- do.call(cbind, coef)
    ## Need to fix this line
    ifelse(Inter,
           rownames(coef) <- c("(Intercept)", colnames(X)),
           rownames <- colnames(X))
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
    df <- do.call(rbind, df)
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
    colnames(df) <- c("model", "variance")
    res <- list(automatic = automatic, call = cl, coef = cbind(drop(coef)), df = df, Inter = Inter, isScaled = isScaled, lambda = lambda, scales = Xscale, terms = Terms,  x = X, xm = Xm, y = Y, model_frame = m)
    ##
    if(!is.null(nPCs))
      {
        if(automatic)
          {
            res$max.nPCs <- max.nPCs
          }
        res$chosen.nPCs <- chosen.nPCs
      }
    class(res) <- "ridgeLogistic"
    res
  }
