rootogram <- function(x, ...)
{
  UseMethod("rootogram")
}

rootogram.goodfit <- function(x, ...)
{
  rootogram.default(x$observed, x$fitted, names = x$count, ...)
}

rootogram.default <- function(x, fitted, names = NULL, scale = c("sqrt", "raw"),
                              type = c("hanging", "standing", "deviation"),
			      bar.col = grey(0.7), line.col = 2,
			      xlab = NULL, ylab = NULL, ylim = NULL, ...)
{
   if(is.null(names)) names <- names(x)
   if(is.table(x)) {
     if(length(dim(x)) > 1) stop ("x must be a 1-way table")
     x <- as.vector(x)
   }
   obs <- x
   fit <- fitted
   if(is.null(xlab)) {xlab <-  "Number of Occurrences"}

   if(match.arg(scale) == "sqrt") {
     obs <- sqrt(obs)
     fit <- sqrt(fit)
     if(is.null(ylab)) {ylab <- "sqrt(Frequency)"}
   } else {
     if(is.null(ylab)) {ylab <- "Frequency"} }


   switch(match.arg(type),

   "hanging" = {
     if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                        c(fit-obs,fit)) + c(0, 0.1)}
     dummy <- barplot(obs, names = names, col = bar.col, beside = FALSE,
             xlab = xlab, ylab = ylab, shift = fit - obs, ylim = ylim, ...)
     lines(dummy, fit, col = line.col, type = "b", pch = 19)
     abline(h = 0)
   },

   "standing" = {
     if(is.null(ylim)) {ylim <- range(-0.01 * c(obs,fit), c(obs,fit)) }
     dummy <- barplot(obs, names = names, col = bar.col,
                      xlab = xlab, ylab = ylab, ylim = ylim, ...)
     lines(dummy, fit, col = line.col, type = "b", pch = 19)
   },

   "deviation" = {
     if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                        c(fit-obs,fit)) + c(0, 0.1)}
     dummy <- barplot(fit - obs, names = names, col = bar.col,
                      xlab = xlab, ylab = ylab, ylim = ylim, ...)
     lines(dummy, fit, col = line.col, type = "b", pch = 19)
   })
}

