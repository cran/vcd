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
           reverse_y = TRUE,
           type = c("observed", "expected"),
           main = deparse(substitute(x)),
           values = c("none", "cells", "margins", "both"),
           frequencies = c("absolute", "relative"),
           sieve_colors = c("red", "blue"),
           sieve_lty = c("longdash", "solid"),
           exp_color = "gray",
           exp_lty = "dotted",
           margin = 0.01,
           newpage = TRUE,
           pop = TRUE,
           margins = c(4,3,4,4),
           xlab = names(dimnames(x))[2],
           ylab = names(dimnames(x))[1],
           ...)
{
  ## parameter handling
  if (length(dim(x)) > 2)
    stop ("Function only implemented for two-way tables")
  ## main
  
  type <- match.arg(type)
  values <- match.arg(values)
  frequencies <- match.arg(frequencies)
  if (is.null(main))
      main <- if (type == "observed") "Sieve diagram" else "Expected frequencies"

  nc <- ncol(x)
  nr <- nrow(x)
  if (reverse_y) x <- x[nr:1,]

  ## compute relative frequencies
  n <- sum(x)
  colFreqs <- colSums(x) / n
  rowFreqs <- rowSums(x) / n

  ## expected values
  ex <- rowFreqs %o% colFreqs * n

  ## signs of deviations
  sgn <- ex - x < 0

  if (newpage) grid.newpage()
  if(!is.null(main))
    margins[2] <- margins[2] + 2
  
  if (values %in% c("margins", "both")) {
    margins[3] <- margins[3] + 1
    margins[4] <- margins[4] + 1
  }
  pushViewport(plotViewport(margins))
  pushViewport(viewport(w = unit(1, "snpc"), h = unit(1, "snpc")))

  if(!is.null(main))
     grid.text(main, y = unit(1.1, "npc"),
               gp = gpar(fontsize = 25))

  
  ## box coordinates for expected areas
  x1 <- c(0, cumsum(colFreqs + margin)[-nc])
  x2 <- x1 + colFreqs
  xmid <- (x1 + x2) / 2
  
  y2 <- 1 + (nr - 1) * margin - c(0, cumsum(rowFreqs + margin)[-nr]) 
  y1 <- y2 - rowFreqs
  ymid <- (y1 + y2) / 2

  ## axis labels
  grid.text(xlab, y = -0.1, gp = gpar(fontsize = 20))
  grid.text(ylab, x = -0.1, gp = gpar(fontsize = 20), rot = 90)
  
  is <- 1:nr
  js <- 1:nc
  
  ## labels
  grid.text(dimnames(x)[[1]][is], x = -0.03, y = ymid[is],
            rot = 90, check.overlap = TRUE, ...)
  
  ## optionally, write marginal frequencies
  if (values %in%  c("margins", "both"))
    grid.text(if (frequencies == "relative") round(rowFreqs[is], 2)
    else round(rowFreqs[is] * n, 1),
              x = 1 + nc * margin + 0.02, y = ymid[is],
              gp = gpar(fontsize = 12, fontface = 2),
              rot = 90, ...)
  
  grid.text(dimnames(x)[[2]][js], x = xmid[js], 
            y = -0.03 + values %in% c("margins", "both") *
            (1 + nr * margin + 0.04), check.overlap = TRUE, ...)
  
  ## optionally, write marginal frequencies
  if (values %in%  c("margins","both"))
    grid.text(if (frequencies == "relative") round(colFreqs[js], 2)
    else round(colFreqs[js] * n, 1),
              x = xmid[js], y = -0.03,
              gp = gpar(fontsize = 12, fontface = 2),  ...)

  ## compute grid
  ltype <- lcolor <- coord <- vector("list", nc * nr)
  cc <- 0
  for (i in 1:nr)
    for (j in 1:nc) {
      cc <- cc + 1
      
      dev <- sgn[i, j] + 1
      line.color <- if (type == "observed") sieve_colors[dev] else exp_color
      line.type  <- if (type == "observed") sieve_lty[dev] else exp_lty
      
      square.side <- sqrt(colFreqs[j] * rowFreqs[i] /
                          ifelse (type == "observed", x[i, j], ex[i, j]))
      ii <- seq(0, rowFreqs[i], by = square.side)
      jj <- seq(0, colFreqs[j], by = square.side)
      
      coord[[cc]] <-
        rbind(cbind(x1[j], x2[j], y1[i] + ii, y1[i] + ii),
              cbind(x1[j] + jj, x1[j] + jj, y1[i], y2[i])
              )
      lcolor[[cc]] <- rep(line.color, length(ii) + length(jj))
      ltype[[cc]] <- rep(line.type, length(ii) + length(jj))
    }

  ## draw grid
  coord <- do.call("rbind", coord)
  ltype <- unlist(ltype)
  lcolor <- unlist(lcolor)
  grid.segments(coord[,1], coord[,3], coord[,2], coord[,4],
                gp = gpar(col = lcolor, lty = ltype),
                default = "native"
                )
  
  ## border
  is <- rep(1:nr, times = nc)
  js <- rep(1:nc, each  = nr)
  grid.rect(x1[js], y1[is], x2[js] - x1[js], y2[is] - y1[is],
            just = c("left", "bottom"), default = "native", ...)
  
  ## optionally, write cell frequencies
  if (values %in% c("cells", "both"))
    grid.text(if (frequencies == "relative")
                round((if (type == "observed") x[is, js] else ex[is, js]) / n, 2)
              else
                round((if (type == "observed") x[is, js] else ex[is, js]), 1),
              xmid[js], ymid[is], gp = gpar(fontsize = 12, fontface = 2),
              check.overlap = TRUE, ...
              )
  if (pop) popViewport(2) else upViewport(2)
  invisible(x)
}





