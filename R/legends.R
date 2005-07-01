legend_resbased <- function(fontsize = 12,
                            x = unit(1, "lines"),
                            y = unit(0.1, "npc"),
                            height = unit(0.8, "npc"),
                            width = unit(0.7, "lines"),
			    digits = 3, 
			    check_overlap = TRUE,
                            text = NULL,
                            steps = 200,
                            ticks = 10,
                            pvalue = TRUE) {

  if(!is.unit(x)) x <- unit(x, "native")
  if(!is.unit(y)) y <- unit(y, "npc")
  if(!is.unit(width)) width <- unit(width, "lines")
  if(!is.unit(height)) height <- unit(height, "npc")
  
  function(residuals, shading, autotext) {
    res <- as.vector(residuals)
    
    if(is.null(text)) text <- autotext
    
    pushViewport(viewport(x = x, y = y, just = c("left", "bottom"),
                          yscale = range(res), default.unit = "native",
                          height = height, width = width))

    p.value <- attr(shading, "p.value")
    legend <- attr(shading, "legend")
    
    if(is.null(legend$col.bins)) {
      col.bins <- seq(min(res), max(res), length = steps)
      at <- NULL
    } else {
      col.bins <- sort(unique(c(legend$col.bins, range(res))))
      col.bins <- col.bins[col.bins <= max(res) & col.bins >= min(res)]
      at <- col.bins
    }
    y.pos <- col.bins[-length(col.bins)]
    y.height <- diff(col.bins)

    grid.rect(x = unit(rep.int(0, length(y.pos)), "npc"),
              y = y.pos,
              height = y.height, default.unit = "native",
              gp = gpar(fill = shading(y.pos + 0.5 * y.height)$fill, col = NULL),
              just = c("left", "bottom"))

    grid.rect()

    if(is.null(at)) at <- seq(from = head(col.bins, 1), to = tail(col.bins, 1), length = ticks)
    grid.text(format(signif(at, digits = digits)),
              x = unit(1, "npc") + unit(0.8, "lines") + unit(1, "strwidth", "-4.44"),
              y = at,
              default.unit = "native", just = c("right", "center"), check.overlap = check_overlap)
    grid.segments(x0 = unit(1, "npc"), x1 = unit(1,"npc") + unit(0.5, "lines"),
                  y0 = at, y1 = at, default.unit = "native")

    popViewport(1)
    grid.text(text, x = x, y = unit(1, "npc") - y + unit(1, "lines"),
              gp = gpar(fontsize = fontsize, lineheight = 0.8),
              just = c("left", "bottom")
              )
    if(!is.null(p.value) && pvalue) {
      grid.text(paste("p-value =\n", format.pval(p.value), sep = ""),
                x = x,
                y = y - unit(1, "lines"),
                gp = gpar(fontsize = fontsize, lineheight = 0.8),
                just = c("left", "top"))
    }
  }
}
class(legend_resbased) <- "panel_generator"

legend_fixed <- function(fontsize = 12,
                         x = unit(1, "lines"),
                         y = unit(0.2, "npc"),
                         height = unit(0.8, "npc"),
                         width = unit(1.5, "lines"),
			 digits = 3, 
			 check_overlap = TRUE,
                         text = NULL) {
  
  if(!is.unit(x)) x <- unit(x, "native")
  if(!is.unit(y)) y <- unit(y, "npc")
  if(!is.unit(width)) width <- unit(width, "lines")
  if(!is.unit(height)) height <- unit(height, "npc")

  function(residuals, shading, autotext) {
    res <- as.vector(residuals)

    if(is.null(text)) text <- autotext
    
    pushViewport(viewport(x = x, y = y, just = c("left", "bottom"),
                          yscale = range(res), default.unit = "native",
                          height = height, width = width))

    p.value <- attr(shading, "p.value")
    legend <- attr(shading, "legend")
    
    if(is.null(legend$col.bins)) {
      col.bins <- seq(min(res), max(res), length = steps)
      at <- NULL
    } else {
      col.bins <- sort(unique(c(legend$col.bins, range(res))))
      col.bins <- col.bins[col.bins <= max(res) & col.bins >= min(res)]
      at <- col.bins
    }
    y.pos <- col.bins[-length(col.bins)]
    y.height <- diff(col.bins)

    grid.rect(x = unit(rep.int(0, length(y.pos)), "npc"),
              y = y.pos,
              height = y.height, default.unit = "native",
              gp = gpar(fill = shading(y.pos + 0.5 * y.height)$fill, col = NULL),
              just = c("left", "bottom"))

    grid.rect()

    ## labeling for fixed intervals as returned in shading.objects
    if(is.null(at)) at <- seq(from = head(col.bins, 1), to = tail(col.bins, 1), length = ticks)
    grid.text(format(signif(at, digits = digits)),
              x = unit(1, "npc") + unit(0.8, "lines") + unit(1, "strwidth", "-4.44"),
              y = at,
              default.unit = "native", just = c("right", "center"), check.overlap = check_overlap)
    grid.segments(x0 = unit(0, "npc"), x1 = unit(1,"npc") + unit(0.5, "lines"),
                  y0 = at, y1 = at, default.unit = "native")

    popViewport(1)
    grid.text(text, x = x + 0.5 * width, y = 0,
              gp = gpar(fontsize = fontsize, lineheight = 0.8),
              just = c("left", "top"),
              rot = 90
              )
  }
}
class(legend_fixed) <- "panel_generator"
