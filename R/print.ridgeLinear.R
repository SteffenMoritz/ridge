## print function for ridgeLinear objects

#' @rdname print
#' @export
#' @importFrom stats coef
print.ridgeLinear <- function(x, all.coef = FALSE, ...)
  {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(coef(x, all.coef = all.coef), ...)
    cat("\n")
    invisible(x)
  }
