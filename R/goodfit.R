goodfit <- function(obj, type = c("poisson", "binomial", "nbinomial"),
                    method = c("ML", "MinChisq"), par = NULL)
{
    if(is.vector(obj)) {
      obj <- table(obj)
    }
    if(is.table(obj)) {
      if(length(dim(obj)) > 1) stop ("obj must be a 1-way table")
      freq <- as.vector(obj)
      count <- as.numeric(names(obj))
    } else {
      if(!(!is.null(ncol(obj)) && ncol(obj) == 2))
        stop("obj must be a 2-column matrix or data.frame")
      freq <- as.vector(obj[,1])
      count <- as.vector(obj[,2])
    }

    ## eliminate zero frequencies
    count <- count[!(freq < 1)]
    freq <- freq[!(freq < 1)]
    n <- length(count)
    df <- n - 1

    type <- match.arg(type)
    method <- match.arg(method)

    switch(type,

    "poisson" = {

      if(!is.null(par)) {
        if(!is.list(par)) stop("`par' must be a named list")
        if(!(names(par) == "lambda")) stop("`par' must specify `lambda'")
	par <- par$lambda
	method <- "fixed"
      }
      else if(method == "ML") {
        df <- df - 1
	par <- weighted.mean(count,freq)
      }
      else if(method == "MinChisq") {
        df <- df - 1

	chi2 <- function(x)
        {
	  p.hat <- diff(c(0, ppois(count[-n], lambda = x), 1))
	  expected <- sum(freq) * p.hat
          sum((freq - expected)^2/expected)
        }

	par <- optimize(chi2, range(count))$minimum
      }
      par <- list(lambda = par)
      p.hat <- dpois(count, lambda = par$lambda)
    },

    "binomial" = {
      size <- par$size
      if(is.null(size)) {
        size <- max(count)
        warning("size was not given, taken as maximum count")
      }

      if(!is.null(par$prob)) {
        if(!is.list(par)) stop("`par' must be a named list and specify `prob'")
	par <- par$prob
	method <- "fixed"
      }
      else if(method == "ML") {
        df <- df - 1
	par <- weighted.mean(count/size, freq)
      }
      else if(method == "MinChisq") {
        df <- df - 1

	chi2 <- function(x)
        {
	  p.hat <- diff(c(0, pbinom(count[-n], prob = x, size = size), 1))
          expected <- sum(freq) * p.hat
          sum((freq - expected)^2/expected)
        }

	par <- optimize(chi2, c(0,1))$minimum
      }
      par <- list(prob = par, size = size)
      p.hat <- dbinom(count, prob = par$prob, size = par$size)
    },


    "nbinomial" = {

      if(!is.null(par)) {
        if(!is.list(par)) stop("`par' must be a named list")
        if(is.character(all.equal(sum(match(names(par), c("size", "prob"))), 3))) stop("`par' must specify `prob' and `size'")
	method <- "fixed"
	par <- c(par$size, par$prob)
      }
      else if(method == "ML") {
        df <- df - 2

	require(package = MASS)
        par <- fitdistr(rep(count, freq), "negative binomial")$estimate
        par <- par[1]/c(1, sum(par))
     }
     else if(method == "MinChisq") {
       df <- df - 2

	## MM
	xbar <- weighted.mean(count,freq)
	s2 <- var(rep(count,freq))
	p <- xbar / s2
	size <- xbar^2/(s2 - xbar)
        par1 <- c(size, p)

	## minChisq
        chi2 <- function(x)
        {
	  p.hat <- diff(c(0, pnbinom(count[-n], size = x[1], prob = x[2]), 1))
          expected <- sum(freq) * p.hat
          sum((freq - expected)^2/expected)
        }

	par <- optim(par1, chi2)$par
      }
      par <- list(size = par[1], prob = par[2])
      p.hat <- dnbinom(count, size = par$size, prob = par$prob)
    })

    expected <- sum(freq) * p.hat

    RVAL <- list(observed = freq,
                 count = count, fitted = expected,
		 type = type, method = method, df = df,
		 par = par)
    class(RVAL) <- "goodfit"
    RVAL
}

print.goodfit <- function(x, ...)
{
    cat(paste("\nObserved and fitted values for", x$type, "distribution\n"))
    if(x$method == "fixed")
      cat("with fixed parameters \n\n")
    else
      cat(paste("with paramaters estimated by `", x$method, "' \n\n", sep = ""))
    RVAL <- cbind(x$count, x$observed, x$fitted)
    colnames(RVAL) <- c("count", "observed", "fitted")
    rownames(RVAL) <- rep("", nrow(RVAL))
    print(RVAL)
    invisible(x)
}

summary.goodfit <- function(object, ...)
{
    df <- object$df
    obsrvd <- object$observed
    count <- object$count
    expctd <- fitted(object)

    G2 <- sum(obsrvd * log(obsrvd/expctd)) * 2

    n <- length(obsrvd)
    switch(object$type,
    "poisson" = { pfun <- "ppois" },
    "binomial" = { pfun <- "pbinom" },
    "nbinomial" = { pfun <- "pnbinom" })
    p.hat <- diff(c(0, do.call(pfun, c(list(q = count[-n]), object$par)), 1))
    expctd <- p.hat * sum(obsrvd)
    X2 <- sum((obsrvd - expctd)^2/expctd)

    names(G2) <- "Likelihood Ratio"
    names(X2) <- "Pearson"
    if(any(expctd) < 5 & object$method != "ML") warning("Chi-squared approximation may be incorrect")

    switch(object$method,
    "ML" = { RVAL <- G2 },
    "MinChisq" = { RVAL <- X2 },
    "fixed" = { RVAL <- c(X2, G2) })

    RVAL <- cbind(RVAL, df, pchisq(RVAL, df = df, lower = FALSE))
    colnames(RVAL) <- c("X^2", "df", "P(> X^2)")

    cat(paste("\n\t Goodness-of-fit test for", object$type, "distribution\n\n"))
    print(RVAL)
    invisible(RVAL)
}

plot.goodfit <- function(x, ...)
{
  rootogram(x, ...)
}

fitted.goodfit <- function(object, ...)
{
  object$fitted
}

predict.goodfit <- function(object, newcount = NULL, type = c("response", "prob"), ...)
{
  if(is.null(newcount)) newcount <- object$count
  type <- match.arg(type)

  switch(object$type,
  "poisson" = { densfun <- "dpois" },
  "binomial" = { densfun <- "dbinom" },
  "nbinomial" = { densfun <- "dnbinom" })

  RVAL <- do.call(densfun, c(list(x = newcount), object$par))
  if(type == "response") RVAL <- RVAL * sum(object$observed)
  return(RVAL)
}

