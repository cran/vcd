assocstats <- function(x) {
  if(!is.matrix(x))
    stop("Function only defined for 2-d - tables.")
  
  tab    <- summary(loglm(~1+2, x))$tests
  phi    <- sqrt(tab[2,1] / sum(x))
  cont   <- sqrt(phi^2 / (1 + phi^2))
  cramer <- sqrt(phi^2 / min(dim(x) - 1))
  structure(
            list(table = x,
                 chisq_tests = tab,
                 phi = phi,
                 contingency = cont,
                 cramer = cramer),
            class = "assocstats"
            )
}

print.assocstats <- function(x,
                             digits = 3,
                             ...)
{
  print(x$chisq_tests, digits = 5, ...)
  cat("\n")
  cat("Phi-Coefficient   :", round(x$phi,    digits = digits), "\n")
  cat("Contingency Coeff.:", round(x$cont,   digits = digits), "\n")
  cat("Cramer's V        :", round(x$cramer, digits = digits), "\n")
  invisible(x)
}

summary.assocstats <- function(object, percentage = FALSE, ...) {
  tab <- summary(object$table, percentage = percentage, ...)
  tab$chisq <- NULL
  structure(list(summary = tab,
                 object  = object),
            class   = "summary.assocstats"
            )
}

print.summary.assocstats <- function(x, ...) {
  cat("\n")
  print(x$summary, ...)
  print(x$object, ...)
  cat("\n")
  invisible(x)
}

