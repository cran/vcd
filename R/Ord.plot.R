Ord.plot <- function(obj, legend = TRUE, estimate = TRUE, tol = 0.1,
                     type = NULL, ylim = NULL, xlab = "Number of occurrences",
		     ylab = "Frequency ratio", main = "Ord plot", ...)
{
  if(is.vector(obj)) {
    obj <- table(obj)
  }
  if(is.table(obj)) {
    if(length(dim(obj)) > 1) stop ("obj must be a 1-way table")
    x <- as.vector(obj)
    count <- as.numeric(names(obj))
  } else {
    if(!(!is.null(ncol(obj)) && ncol(obj) == 2))
      stop("obj must be a 2-column matrix or data.frame")
    x <- as.vector(obj[,1])
    count <- as.vector(obj[,2])
  }

  y <- count * x/c(NA, x[-length(x)])
  fm <- lm(y ~ count)
  fmw <- lm(y ~ count, weights = sqrt(x - 1))
  fit1 <- predict(fm, data.frame(count))
  fit2 <- predict(fmw, data.frame(count))
  if(is.null(ylim)) ylim <- range(c(y, fit1, fit2), na.rm = TRUE)
  plot(y ~ count, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  lines(count, fit1)
  lines(count, fit2, col = 2)
  RVAL <- coef(fmw)
  names(RVAL) <- c("Intercept", "Slope")
  if(legend)
  {
    legend.text <- c(paste("slope =", round(RVAL[2], digits = 3)),
                     paste("intercept =", round(RVAL[1], digits = 3)))
    if(estimate) {
      ordfit <- Ord.estimate(RVAL, type = type, tol = tol)
      legend.text <- c(legend.text, "", paste("type:", ordfit$type),
        paste("estimate:", names(ordfit$estimate),"=", round(ordfit$estimate, digits = 3)))
    }
    legend(min(count), ylim[2], legend.text, bty = "n")
  }
  invisible(RVAL)
}

Ord.estimate <- function(x, type = NULL, tol = 0.1)
{
  a <- x[1]
  b <- x[2]
  if(!is.null(type))
    type <- match.arg(type, c("poisson", "binomial", "nbinomial", "log-series"))
  else {
    if(abs(b) < tol) type <- "poisson"
    else if(b < (-1 * tol)) type <- "binomial"
    else if(a > (-1 * tol)) type <- "nbinomial"
    else if(abs(a + b) < 4*tol) type <- "log-series"
    else type <- "none"
  }

  switch(type,

  "poisson" = {
    par <- a
    names(par) <- "lambda"
    if(par < 0) warning("lambda not > 0")
  },
  "binomial" = {
    par <- b/(b - 1)
    names(par) <- "prob"
    if(abs(par - 0.5) > 0.5) warning("prob not in (0,1)")
  },
  "nbinomial" = {
    par <- 1 - b
    names(par) <- "prob"
    if(abs(par - 0.5) > 0.5) warning("prob not in (0,1)")
  },

  "log-series" = {
    par <- b
    names(par) <- "theta"
    if(par < 0) warning("theta not > 0")
  },
  "none" = {
    par <- NA
  })
  list(estimate = par, type = type)
}


