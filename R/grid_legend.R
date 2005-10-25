grid_legend <-
function (x, y, pch, col, labels, frame = TRUE, hgap = unit(0.5, "lines"), 
    vgap = unit(0.3, "lines"), default_units = "lines", gp = gpar(), 
    draw = TRUE, title = "Legend:") 
{
    labels <- as.character(labels)
    if (is.logical(title) && !title)
      title <- NULL
    if (!is.null(title)) {
      labels <- c(title, labels)
      pch <- c(NA, pch)
      col <- c(NA, col)
    }
    nkeys <- length(labels)
    if (length(pch) != nkeys) 
        stop("pch and labels not the same length")
    if (!is.unit(hgap)) 
        hgap <- unit(hgap, default_units)
    if (length(hgap) != 1) 
        stop("hgap must be single unit")
    if (!is.unit(vgap)) 
        vgap <- unit(vgap, default_units)
    if (length(vgap) != 1) 
        stop("vgap must be single unit")
    legend.layout <- grid.layout(nkeys, 3, widths = unit.c(unit(2, 
        "lines"), max(unit(rep(1, nkeys), "strwidth", as.list(labels))), 
        hgap), heights = unit.pmax(unit(1, "lines"), vgap + unit(rep(1, 
        nkeys), "strheight", as.list(labels))))
    fg <- frameGrob(layout = legend.layout, gp = gp)
    for (i in 1:nkeys) {
      tit <- !is.null(title) && i == 1
      if (!tit)
        fg <- placeGrob(fg, pointsGrob(0.5, 0.5, pch = pch[i], gp = gpar(col = col[i])), 
                        col = 1, row = i)
      fg <- placeGrob(fg, textGrob(labels[i],
                                   x = 0 + 0.3 * tit,
                                   y = 0.5, 
                                   just = c("left", "center")),
                      col = 2 - tit, row = i)
    }
    pushViewport(viewport(x, y, height = unit(nkeys, "lines"),
                          width = grobWidth(fg)))
    if (draw) 
      grid.draw(fg)
    if (frame)
      grid.rect(gp = gpar(fill = "transparent"))
    popViewport(1)
    invisible(fg)
}
