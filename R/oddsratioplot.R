"oddsratio" <-
function (x, stratum = NULL, log = TRUE, conf.level = 0.95) {
  l <- length(dim(x))
  if (l > 2 && is.null(stratum))
    stratum <- 3:l
  if (l - length(stratum) > 2)
    stop("All but 2 dimensions must be specified as strata.")
  if (l == 2 && dim(x) != c(2,2))
    stop("Not a 2 x 2 - table.")
  if (!is.null(stratum) && dim(x)[-stratum] != c(2,2))
    stop("Need strata of 2 x 2 - tables.")
 
  lor <- function (y) {
    y <- y + 0.5 
    or <- y[1,1] * y[2,2] / y[1,2] / y[2,1]
    if (log) log(or) else or
  }

  ase <- function(y) sqrt(sum(1/(y + 0.5)))

  if(is.null(stratum)) {
    LOR <- lor(x)
    ASE <- ase(x)
  } else {
    LOR <- apply(x, stratum, lor)
    ASE <- apply(x, stratum, ase)
  }

  I <- ASE * qnorm((1 + conf.level) / 2)
  Z <- LOR / ASE
  
  structure(LOR,
            ASE = if(log) ASE,
            lwr = if(log) LOR - I else exp(log(LOR) - I),
            upr = if(log) LOR + I else exp(log(LOR) + I),
            Z   = if(log) Z,
            P   = if(log) 1 - pnorm(abs(Z)),
            log = log,
            class = "oddsratio"
            )}

"print.oddsratio" <-
function(x, ...) {
  if (length(dim(x)) > 1)
    print(cbind(unclass(x)))
  else
    print(c(x))
  invisible(x)
}

"summary.oddsratio" <-
function(object, ...) {
  if(!is.null(dim(object)))
    ret <- object
  else {
    ret <- cbind(object,
          ASE = attr(object, "ASE"),
          Z   = attr(object, "Z"),
          P   = attr(object, "P"),
          lwr = attr(object, "lwr"),
          upr = attr(object, "upr")
          )
    colnames(ret)[1] <- if(attr(object, "log")) "Log Odds Ratio" else "Odds Ratio"
  }
  
  class(ret) <- "summary.oddsratio"
  ret
}


"print.summary.oddsratio" <-
function(x, ...) {
  if(!is.null(attr(x, "log"))) {
    cat("\n")
    cat(if(attr(x, "log")) "Log Odds Ratio(s):" else "Odds Ratio(s):", "\n\n")
    print.oddsratio(x)
    cat("\nAsymptotic Standard Error(s):\n\n")
    print(attr(x, "ASE"))
    cat("\n")
  } else print(unclass(x))
  invisible(x)
}

"plot.oddsratio" <-
function(x,
         confidence = TRUE,
         type = "o",
         ylab = NULL,
         xlab = "Strata",
         whiskers = 0.1,
         ...)
{
  if (length(dim(x)) > 1)
    stop ("Plot function works only on vectors.")
  
  yrange <- range(x)
  
  if(confidence) {
    lwr <- attr(x, "lwr")
    upr <- attr(x, "upr")
    yrange[1] <- trunc(min(yrange[1], min(lwr)))
    yrange[2] <- ceiling(max(yrange[2], max(upr)))
  }

  plot(unclass(x),
       xlab = xlab,
       ylab = if(!is.null(ylab)) ylab else if(attr(x, "log")) "Log Odds Ratio" else "Odds Ratio",
       type = type,
       xaxt = "n",
       ylim = yrange,
       ...)
  axis (1, at = 1:length(x), names(x))

  if (confidence)
    for (i in 1:length(x)) {
      lines(c(i, i), c(lwr[i], upr[i]))
      lines(c(i - whiskers/2, i + whiskers/2), c(lwr[i], lwr[i]))
      lines(c(i - whiskers/2, i + whiskers/2), c(upr[i], upr[i]))
    }
}









