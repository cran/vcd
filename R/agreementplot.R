"agreementplot" <- function (x, ...)
  UseMethod ("agreementplot")

"agreementplot.formula" <-
function (formula, data = NULL, ..., subset) 
{
    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())
    if (inherits(edata, "ftable") || inherits(edata, "table")) {
        data <- as.table(data)
        varnames <- attr(terms(formula), "term.labels")
        if (all(varnames != ".")) 
            data <- margin.table(data, match(varnames, names(dimnames(data))))
        agreementplot(data, ...)
    }
    else {
        if (is.matrix(edata)) 
            m$data <- as.data.frame(data)
        m$... <- NULL
        m[[1]] <- as.name("model.frame")
        mf <- eval(m, parent.frame())
        if (length(formula) == 2) {
          by <- mf
          y <- NULL
        }
        else {
          i <- attr(attr(mf, "terms"), "response")
          by <- mf[-i]
          y <- mf[[i]]
        }
        by <- lapply(by, factor)
        x <- if (is.null(y)) 
          do.call("table", by)
        else if (NCOL(y) == 1) 
          tapply(y, by, sum)
        else {
          z <- lapply(as.data.frame(y), tapply, by, sum)
          array(unlist(z), dim = c(dim(z[[1]]), length(z)), dimnames = c(dimnames(z[[1]]), 
                                                              list(names(z))))
        }
        x[is.na(x)] <- 0
        agreementplot(x, ...)
    }
}

"agreementplot.default" <-
  function(x,
           reverse.y = TRUE,
           main = "Agreement Chart",
           weights = c(1, 1 - 1 / (ncol(x) - 1)^2),
           cex.main = 2,
           cex.lab = 1.5,
           xlab = names(dimnames(x))[2],
           ylab = names(dimnames(x))[1],
           ...)
{
  if (length(dim(x)) > 2)
    stop("Function implemented for two-way tables only!")
  if (ncol(x) != nrow(x))
    stop("Dimensions must have equal length!")

  nc <- ncol(x)
  
  ## compute relative frequencies
  n <- sum(x)
  colFreqs <- colSums(x) / n
  rowFreqs <- rowSums(x) / n

  ## margins, limits (hard-coded, argh!)
  bm <- 0.1
  lm <- 0.1
  tm <- 0.1
  rm <- 0.1

  xlim = c(0, 1 + lm + rm)
  ylim = c(0, 1 + tm + bm)

  ## init device
  opar <- par(usr = c(xlim, ylim), mar = c(0, 0, 0, 0))
  on.exit(par(opar))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1)

  ## title
  text(x = lm + 1 / 2, y = ylim[2], labels = main, cex = cex.main)

  ## axis labels
  text(x = lm + 1 / 2, y = 0, labels = xlab, cex = cex.lab)
  text(x = 0, y = bm + 1 / 2, labels = ylab, cex = cex.lab, srt = 90)
  
  rect(lm, bm, lm + 1, bm + 1)

  xc <- c(0, cumsum(colFreqs))
  yc <- c(0, cumsum(rowFreqs))

  my.text <- function (y, ...)
    if (reverse.y)
      text(y, ...)
    else
      text(2 * bm + 1 - y, ...)

  my.rect <- function (xleft, ybottom, xright, ytop, ...)
    if (reverse.y)
      rect(lm + xleft, bm + ybottom, lm + xright, bm + ytop, ...)
    else
      rect(lm + xleft, 1 + tm - ybottom, lm + xright, 1 + tm - ytop, ...)
  
  A <- matrix(0, length(weights), nc)
  for (i in 1:nc) {
    ## x - axis
    text(x = lm + xc[i] + (xc[i+1] - xc[i]) / 2, y = bm - 0.04,
         labels = dimnames(x)[[2]][i], ...)

    ## y - axis
    my.text(y = bm + yc[i] + (yc[i+1] - yc[i]) / 2, x = lm - 0.03,
            labels = dimnames(x)[[1]][i], srt = 90, ...)
    
    ## expected rectangle
    my.rect(xc[i], yc[i], xc[i+1], yc[i+1])
    
    ## observed rectangle
    y0 <- c(0, cumsum(x[i,])) / sum(x[i,])
    x0 <- c(0, cumsum(x[,i])) / sum(x[,i])

    rec <- function (col, dens)
      my.rect(xc[i] + (xc[i+1] - xc[i]) * x0[lb],
              yc[i] + (yc[i+1] - yc[i]) * y0[lb],
              xc[i] + (xc[i+1] - xc[i]) * x0[tr],
              yc[i] + (yc[i+1] - yc[i]) * y0[tr],
#             col = gray(1-(weights[j])^2)
              col = col,
              density = dens,
              angle = 135
              )

    for (j in length(weights):1) {
      lb <- max(1, i - j + 1)
      tr <- 1 + min(nc, i + j - 1)
      A[j, i] <- sum(x[lb:(tr-1),i]) * sum(x[i, lb:(tr-1)])
      rec("white", NULL) ## erase background
      rec("black", if (weights[j] < 1) weights[j] * 20 else NULL)
    }

    ## correct A[j,i] -> not done by Friendly==Bug?
    for (j in length(weights):1) 
      if (j > 1) A[j, i] <- A[j, i] - A[j - 1, i]
  }
  if (reverse.y)
    lines(c(lm, bm + 1), c(lm, 1 + bm), col = "red", lty = "longdash")
  else
    lines(c(lm, bm + 1), c(lm + 1, bm), col = "red", lty = "longdash")
  
  ## Statistics - Returned invisibly
  ads <- crossprod(diag(x)) 
  ar  <- n * n * crossprod(colFreqs, rowFreqs)
  invisible(list(
                 Bangdiwala = ads / ar,
                 Bangdiwala.Weighted = (sum(weights * A)) /  ar,
                 weights = weights,
                 )
            )
}

Kappa <- function (x, weights = c("Equal-Spacing", "Fleiss-Cohen"), conf.level = 0.95)
{
  if (is.character(weights))
      weights = match.arg(weights)

  q <- qnorm((1 + conf.level) / 2)
  
  d  <- diag(x)
  n  <- sum(x)
  nc <- ncol(x)
  colFreqs <- colSums(x)/n
  rowFreqs <- rowSums(x)/n
  
  ## Kappa
  kappa <- function (po, pc)
    (po - pc) / (1 - pc)
  std  <- function (po, pc, W = 1)
    sqrt(sum(W * W * po * (1 - po)) / crossprod(1 - pc) / n)
    
  ## unweighted
  po <- sum(d) / n
  pc <- crossprod(colFreqs, rowFreqs)
  k <- kappa(po, pc)
  s <- std(po, pc)
  
  ## weighted 
  W <- if (is.matrix(weights))
    weights
  else if (weights == "Equal-Spacing")
    outer (1:nc, 1:nc, function(x, y) 1 - abs(x - y) / (nc - 1))
  else
    outer (1:nc, 1:nc, function(x, y) 1 - (abs(x - y) / (nc - 1))^2)
  pow <- sum(W * x) / n
  pcw <- sum(W * colFreqs %o% rowFreqs)
  kw <- kappa(pow, pcw)
  sw <- std(x / n, 1 - pcw, W)

  structure(
            list(Unweighted = c(
                   value = k,
                   ASE   = s,
                   lwr   = k - s * q,
                   upr   = k + s * q
                   ),
                 Weighted = c(
                   value = kw,
                   ASE   = sw,
                   lwr   = kw - sw * q,
                   upr   = kw + sw * q 
                   ),
                 Weights = W
                 ),
            class = "Kappa"
       )
}

print.Kappa <- function (x, ...) {
  tab <- rbind(x$Unweighted, x$Weighted)
  rownames(tab) <- names(x)[1:2]
  print(tab)
  invisible(x)
}

summary.Kappa <- function (object, ...)
  structure(object, class = "summary.Kappa")

print.summary.Kappa <- function (x, ...) {
  print.Kappa(x)
  cat("\nWeights:\n")
  print(x$Weights)
  invisible(x)
}

expected <- function(x, frequency = c("absolute","relative")) {
  if (!is.array(x))
    stop("Need array of absolute frequencies!")
  frequency <- match.arg(frequency)

  n <- sum(x)
  x <- x / n
  d <- length(dim(x))
  tab <- apply(x, 1, sum)
  for (i in 2:d)
    tab <- tab %o% apply(x, i, sum)
  if (frequency == "relative") tab else tab * n
}

mar.table <- function(x) {
  if(!is.matrix(x))
    stop("Function only defined for m x n - tables.")
  tab <- rbind(cbind(x, TOTAL = rowSums(x)), TOTAL = c(colSums(x), sum(x)))
  names(dimnames(tab)) <- names(dimnames(x))
  tab
}

summary.table <- function(object,
                    margins = TRUE,
                    percentages = FALSE,
                    conditionals = c("none", "row", "column"),
                    ...
                    )
{
  ret <- list()
  ret$chisq <- base::summary.table(object, ...)
  
  if(is.matrix(object)) {
    
    conditionals <- match.arg(conditionals)
  
    tab <- array(0, c(dim(object) + margins, 1 + percentages + (conditionals != "none")))

    ## frequencies
    tab[,,1] <- if(margins) mar.table(object) else object

    ## percentages
    if(percentages) {
      tmp <- prop.table(object)
      tab[,,2] <- 100 * if(margins) mar.table(tmp) else tmp
    }

    ## conditional distributions
    if(conditionals != "none") {
      tmp <- prop.table(object, margin = 1 + (conditionals == "column"))
      tab[,,2 + percentages] <- 100 * if(margins) mar.table(tmp) else tmp
    }

    ## dimnames
    dimnames(tab) <- c(dimnames(if(margins) mar.table(object) else object),
                       list(c("freq",
                              if(percentages) "%",
                              switch(conditionals, row = "row%", column = "col%")
                              )
                            )
                       )

    ## patch row% / col% margins
    if(conditionals == "row") 
      tab[dim(tab)[1],,2 + percentages] <- NA
    
    if(conditionals == "column")
      tab[,dim(tab)[2],2 + percentages] <- NA
    
    ret$table <- tab
  }    

  class(ret) <- "summary.table"
  ret
}

print.summary.table <- 
function (x, digits = max(1, getOption("digits") - 3), ...) 
{
  if (!is.null(x$table))
    if(dim(x$table)[3] == 1)
      print(x$table[,,1], digits = digits)
    else
      print(ftable(aperm(x$table, c(1,3,2))), 2, digits = digits)
  
  cat("\n")
  
  if (!is.null(x$chisq))
    base::print.summary.table(x$chisq, digits, ...)
  invisible(x)
}

assoc.stats <- function(x) {
  if(!is.matrix(x))
    stop("Function only defined for m x n - tables.")
  require(MASS)
  
  tab    <- summary(loglm(~1+2, x))$tests
  phi    <- sqrt(tab[2,1] / sum(x))
  cont   <- sqrt(phi^2 / (1 + phi^2))
  cramer <- sqrt(phi^2 / min(dim(x) - 1))
  structure(
            list(table = x,
                 chisq.tests = tab,
                 phi = phi,
                 contingency = cont,
                 cramer = cramer),
            class = "assoc.stats"
            )
}

print.assoc.stats <- function(x,
                              digits = 3,
                              ...)
{
  print(x$chisq.tests, digits = 5)
  cat("\n")
  cat("Phi-Coefficient   :", round(x$phi,    digits = digits), "\n")
  cat("Contingency Coeff.:", round(x$cont,   digits = digits), "\n")
  cat("Cramer's V        :", round(x$cramer, digits = digits), "\n")
  invisible(x)
}

summary.assoc.stats <- function(object, percentage = FALSE, ...) {
  tab <- summary(object$table, percentage = percentage, ...)
  tab$chisq <- NULL
  structure(list(summary = tab,
                 object  = object),
            class   = "summary.assoc.stats"
            )
}

print.summary.assoc.stats <- function(x, ...) {
  cat("\n")
  print(x$summary)
  print(x$object)
  cat("\n")
  invisible(x)
}

woolf.test <- function(x) {
  DNAME <- deparse(substitute(x))
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  o <- log(or)
  e <- weighted.mean(log(or), w)
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Woolf-test on Homogeneity of Odds Ratios (no 3-Way assoc.)"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = o, 
                 expected = e), class = "htest")
}




