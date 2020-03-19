## print function for ridgeLinear objects

#' @rdname print
#' @export
print.ridgeLinear <- function(x, all.coef = FALSE, ...)
  {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(stats::coef(x, all.coef = all.coef), ...)
    cat("\n")
    invisible(x)
  }
