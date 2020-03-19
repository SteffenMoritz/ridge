## S3 method for class 'summary.ridgeLinear'

#' @rdname print
#' @export
#' @importFrom stats printCoefmat
print.summary.ridgeLinear <- function(x, digits = max(3, getOption("digits") - 3),
                                      signif.stars = getOption("show.signif.stars"), ...)
  {
    summaries <- x$summaries
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if(x$automatic && !x$all.coef)
      {
        chosenSummary <- summaries[[x$chosen.nPCs]]
        cat("\nCoefficients:\n")
        coefs <- chosenSummary$coefficients
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                     na.print = "NA", ...)
        cat("\nRidge parameter:", chosenSummary$lambda)
        if(x$automatic)
          cat(", chosen automatically")
        if(!is.null(x$chosen.nPCs))
          cat(gettextf(", computed using %d PCs\n", x$chosen.nPCs))
        else
          cat("\n")
        ## df: degrees of freedom: model, variance, residual
        cat("\nDegrees of freedom: model", format(signif(chosenSummary$df[1], digits)), ", variance", format(signif(chosenSummary$df[2], digits)), ", residual", format(signif(chosenSummary$df[3], digits)), "\n")
        cat("\n")
      } else {
        ## Want to mark out the best chosen lambda
        ## in this bit
        for(i in seq(length(x$summaries)))
          {
            chosenSummary <- summaries[[i]]
            cat("\nCoefficients:\n")
            coefs <- chosenSummary$coefficients
            printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                         na.print = "NA", ...)
            cat("\nRidge parameter:", chosenSummary$lambda)
            if(x$automatic && (x$chosen.nPCs == i))
              cat(", chosen automatically")
            if(!is.null(chosenSummary$nPCs))
              cat(gettextf(", computed using %d PCs\n", chosenSummary$nPCs))
            else
              cat("\n")
            ## df: degrees of freedom: model, variance, residual
            cat("\nDegrees of freedom: model", format(signif(chosenSummary$df[1], digits)), ", variance", format(signif(chosenSummary$df[2], digits)), ", residual", format(signif(chosenSummary$df[3], digits)), "\n")
            cat("\n")
            invisible(x)
          } ## Ends for i in seq(length(x$summaries))
      } ## Ends else 
  }



