distplot <- function(obj, type = c("poisson", "binomial", "nbinomial"),
                     size = NULL, lambda = NULL, legend = TRUE, ylim = NULL,
                     line.col = 2, conf.int = TRUE, conf.level = 0.95, main = NULL,
		     xlab = "Number of occurrences", ylab = "Distribution metameter", ...)
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

  myindex <- (1:length(freq))[freq > 0]
  mycount <- count[myindex]
  myfreq <- freq[myindex]

  switch(match.arg(type),

  "poisson" = {
    par.ml <- goodfit(obj, type = type)$par$lambda

    phi <- function(nk, k, N, size = NULL)
      ifelse(nk > 0, lgamma(k + 1) + log(nk/N), NA)
    y <- phi(myfreq, mycount, sum(freq))
    if(!is.null(lambda)) y <- y + lambda - mycount * log(lambda)
    fm <- lm(y ~ mycount)
    par.estim <- exp(coef(fm)[2])
    names(par.estim) <- "lambda"
    if(!is.null(lambda)) par.estim <- par.estim * lambda
    legend.text <- paste("exp(slope) =", round(par.estim, digits = 3))
    if(is.null(main)) main <- "Poissoness plot"
  },

  "binomial" = {
    if(is.null(size)) {
      size <- max(count)
      warning("size was not given, taken as maximum count")
    }
    par.ml <- goodfit(obj, type = type, par = list(size = size))$par$prob

    phi <- function(nk, k, N, size)
      log(nk) - log(N * choose(size, k))
    y <- phi(myfreq, mycount, sum(freq), size = size)
    fm <- lm(y ~ mycount)
    par.estim <- exp(coef(fm)[2])
    par.estim <- par.estim / (1 + par.estim)
    names(par.estim) <- "prob"
    legend.text <- paste("inv.logit(slope) =", round(par.estim, digits = 3))
    if(is.null(main)) main <- "Binomialness plot"
  },

  "nbinomial" = {
    par.ml <- goodfit(obj, type = type)$par
    size <- par.ml$size
    par.ml <- par.ml$prob
    phi <- function(nk, k, N, size)
      log(nk) - log(N * choose(size + k - 1, k))
    y <- phi(myfreq, mycount, sum(freq), size = size)
    fm <- lm(y ~ mycount)
    par.estim <- 1 - exp(coef(fm)[2])
    names(par.estim) <- "prob"
    legend.text <- paste("1-exp(slope) =", round(par.estim, digits = 3))
    if(is.null(main)) main <- "Negative binomialness plot"
  })

  yhat <- ifelse(myfreq > 1.5, myfreq - 0.67, 1/exp(1))
  yhat <- phi(yhat, mycount, sum(freq), size = size)
  if(!is.null(lambda)) yhat <- yhat + lambda - mycount * log(lambda)

  phat <- myfreq / sum(myfreq)
  ci.width <- qnorm(1-(1 - conf.level)/2) *
              sqrt(1-phat)/sqrt(myfreq - (0.25 * phat + 0.47)*sqrt(myfreq))

  RVAL <- cbind(count, freq, NA, NA, NA, NA, NA)
  RVAL[myindex,3:7] <- cbind(y,yhat,ci.width, yhat-ci.width, yhat + ci.width)
  RVAL <- as.data.frame(RVAL)
  names(RVAL) <- c("Counts", "Freq", "Metameter", "CI.center",
                   "CI.width", "CI.lower", "CI.upper")

  if(is.null(ylim)) ylim <- range(RVAL[,c(3,6,7)], na.rm = TRUE)
  plot(Metameter ~ Counts, ylim = ylim, data = RVAL,
       xlab = xlab, ylab = ylab, main = main, ...)
  abline(fm, col = line.col)

  if(conf.int) {
    points(CI.center ~ Counts, data = RVAL, pch = 19, cex = 0.6)
    arrows(RVAL[,1], RVAL[,6], RVAL[,1], RVAL[,7], length = 0, lty = 3)
  }

  if(legend) {
    mymin <- which.min(RVAL[,5])
    leg.x <- RVAL[mymin,1]
    if(RVAL[mymin,6] - ylim[1] > ylim[2] - RVAL[mymin,7])
      leg.y <- ylim[1] + 0.7 * (RVAL[mymin,6] - ylim[1])
    else leg.y <- ylim[2]

    legend.text <- c(paste("slope =", round(coef(fm)[2], digits = 3)),
                     paste("intercept =", round(coef(fm)[1], digits = 3)),
		     "", paste(names(par.estim),": ML =", round(par.ml, digits=3)),
		     legend.text)
    legend(leg.x, leg.y, legend.text, bty = "n")
  }
  invisible(RVAL)
}
