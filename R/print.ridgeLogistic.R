## print function for ridgeLogistic objects

print.ridgeLogistic <- function(x, all.coef = FALSE, ...)
  {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(coef(x, all.coef = all.coef), ...)
    cat("\n")
    invisible(x)
  }
