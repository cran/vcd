"sieveplot" <- function (x, ...)
  UseMethod ("sieveplot")

"sieveplot.formula" <-
function (formula, data = NULL, ..., subset) 
{
    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())
    if (inherits(edata, "ftable") || inherits(edata, "table")) {
        data <- as.table(data)
        varnames <- attr(terms(formula), "term.labels")
        if (all(varnames != ".")) 
            data <- margin.table(data, match(varnames, names(dimnames(data))))
        sieveplot(data, ...)
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

        sieveplot(x, ...)
      }
}

"sieveplot.default" <-
  function(x,
           reverse.y = TRUE,
           type = c("sieve","expected"),
           main = NULL,
           values = c("none", "cells", "margins", "both"),
           frequencies = c("absolute", "relative"),
           sieve.colors = c("red","blue"),
           sieve.lty = c("longdash", "solid"),
           exp.color = "gray",
           exp.lty = "dotted",
           margin = 0.01,
           cex.main = 2,
           cex.lab = 1.5,
           xlab = names(dimnames(x))[2],
           ylab = names(dimnames(x))[1],
           ...)
{
  ## parameter handling
  if (length(dim(x)) > 2)
    stop ("Function only implemented for two-way tables")
  
  type <- match.arg(type)
  values <- match.arg(values)
  frequencies <- match.arg(frequencies)
  if (is.null(main))
      main <- if (type == "sieve") "Sieve diagram" else "Expected frequencies"

  nc   <- ncol(x)
  nr   <- nrow(x)
  if (reverse.y) x <- x[nr:1,]

  ## compute relative frequencies
  n <- sum(x)
  colFreqs <- colSums(x) / n
  rowFreqs <- rowSums(x) / n

  ## expected values
  ex <- rowFreqs %o% colFreqs * n

  ## signs of deviations
  sgn <- ex - x < 0

  # margins, limits (hard-coded, argh!)
  bm <- 0.1
  lm <- 0.1
  tm <- 0.1 + 0.05 * values %in%  c("margins", "both")
  rm <- 0.1

  xlim <- c(0, 1 + (nc - 1) * margin + lm + rm)
  ylim <- c(0, 1 + (nr - 1) * margin + tm + bm)

  ## box coordinates for expected areas
  x1 <- lm + c(0, cumsum(colFreqs + margin)[-nc])
  x2 <- x1 + colFreqs
  xmid <- (x1 + x2) / 2
  
  y2 <- bm + 1 + (nr - 1) * margin - c(0, cumsum(rowFreqs + margin)[-nr]) 
  y1 <- y2 - rowFreqs
  ymid <- (y1 + y2) / 2

  ## setup device
  opar <- par(mar = c(0, 0, 0, 0))
  on.exit(par(opar))
  plot.new()
  par(usr = c(xlim, ylim))
  plot.window(xlim = xlim, ylim = ylim, asp = 1)
  
  ## title
  text(x = lm + (1 + (nc - 1) * margin) / 2, y = ylim[2], labels = main, cex = cex.main)

  ## axis labels
  text(x = lm + (1 + (nc - 1) * margin) / 2, y = 0, labels = xlab, cex = cex.lab)
  text(x = 0, y = bm + (1 + (nr - 1) * margin) / 2, labels = ylab, cex = cex.lab, srt = 90)
  
  ## boxes
  for (i in 1:nr)
    for (j in 1:nc) {
      
      ## WRITE LABELS

      if (j == 1) { ## y-axis
        text(x = lm - 0.03, y = ymid[i], labels = dimnames(x)[[1]][i], srt = 90, ...)
        ## optionally, write marginal frequencies
        if (values %in%  c("margins", "both"))
          text(x = 1 + nc * margin + lm + 0.03, y = ymid[i], font = 2,
               labels = if (frequencies == "relative") round(rowFreqs[i], 2)
                        else round(rowFreqs[i] * n, 1), srt = 90, ...)
      }
      if (i == 1) { ## x-axis
        text(y = bm - 0.03 + values %in% c("margins", "both") * (1 + nr * margin + 0.04),
             x = xmid[j], labels = dimnames(x)[[2]][j], ...)
        ## optionally, write marginal frequencies
        if (values %in%  c("margins","both"))
          text(pos = 1, x = xmid[j], y = bm - 0.01, font = 2,
               labels = if (frequencies == "relative") round(colFreqs[j], 2)
                        else round(colFreqs[j] * n, 1), ...)
      }

      ## DRAW GRID

      square.side <- sqrt(colFreqs[j] * rowFreqs[i] / if (type == "sieve") x[i, j] else ex[i, j])
      dev <- sgn[i, j] + 1
      line.color <- if (type == "sieve") sieve.colors[dev] else exp.color
      line.type  <- if (type == "sieve") sieve.lty[dev] else exp.lty
      for (ii in seq(0, rowFreqs[i], by = square.side))
        lines(c(x1[j], x2[j]), c(y1[i], y1[i]) + ii,
              col = line.color, lty = line.type
              )
      for (jj in seq(0, colFreqs[j], by = square.side))
        lines(c(x1[j], x1[j]) + jj, c(y1[i], y2[i]),
              col = line.color, lty = line.type
              )
        
      ## OPTIONALLY, WRITE CELL FREQUENCIES

      if (values %in% c("cells", "both"))
          text(xmid[j], ymid[i],
               if (frequencies == "relative")
                 round((if (type == "sieve") x[i, j] else ex[i, j]) / n, 2)
               else
                 round((if (type == "sieve") x[i, j] else ex[i, j]), 1),
               font = 2, ...
               )

      ## BORDER

      rect(x1[j], y1[i], x2[j], y2[i], ...)
    }
}




