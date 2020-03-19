## print function for ridgeLogistic objects

#' @rdname print
#' @export
print.ridgeLogistic <- function(x, all.coef = FALSE, ...)
  {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(stats::coef(x, all.coef = all.coef), ...)
    cat("\n")
    invisible(x)
  }
