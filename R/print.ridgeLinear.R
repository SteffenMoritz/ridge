## print function for ridgeLinear objects

print.ridgeLinear <- function(x, all.coef = FALSE, ...)
  {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(coef(x, all.coef = all.coef), ...)
    cat("\n")
    invisible(x)
  }
