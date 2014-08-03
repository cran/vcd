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
  rect_gp = gpar(fill = "lightgray"), lines_gp = gpar(col = "red"),
  points_gp = gpar(col = "red"), pch = 19,
  xlab = NULL, ylab = NULL, ylim = NULL,
  name = "rootogram", newpage = TRUE, pop = TRUE, ...)
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
     if(is.null(ylab)) {ylab <- "Frequency"}
   }


   switch(match.arg(type),

   "hanging" = {
     if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                        c(fit-obs,fit)) + c(0, 0.1)}
     dummy <- grid_barplot(obs, names = names, offset = fit - obs, gp = rect_gp,
             xlab = xlab, ylab = ylab, ylim = ylim,
	     name = name, newpage = newpage, pop = FALSE, ...)
     downViewport(name)
     grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
     grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
     grid.lines(x = unit(c(0, 1), "npc"), y = unit(0, "native"))
     if(pop) popViewport() else upViewport()
   },

   "standing" = {
     if(is.null(ylim)) {ylim <- range(-0.01 * c(obs,fit), c(obs,fit)) }
     dummy <- grid_barplot(obs, names = names, gp = rect_gp,
             xlab = xlab, ylab = ylab, ylim = ylim,
	     name = name, newpage = newpage, pop = FALSE, ...)
     downViewport(name)
     grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
     grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
     if(pop) popViewport() else upViewport()
   },

   "deviation" = {
     if(is.null(ylim)) {ylim <- range(-0.1 * c(fit-obs,fit),
                        c(fit-obs,fit)) + c(0, 0.1)}
     dummy <- grid_barplot(fit - obs, names = names, gp = rect_gp,
             xlab = xlab, ylab = ylab, ylim = ylim,
	     name = name, newpage = newpage, pop = FALSE, ...)
     downViewport(name)
     grid.lines(x = dummy, y = fit, default.units = "native", gp = lines_gp)
     grid.points(x = dummy, y = fit, default.units = "native", gp = points_gp, pch = pch)
     if(pop) popViewport() else upViewport()
   })
}

grid_barplot <- function(height, width = 0.8, offset = 0,
  names = NULL, xlim = NULL, ylim = NULL, xlab = "", ylab = "", main = "",
  gp = gpar(fill = "lightgray"), name = "grid_barplot", newpage = TRUE, pop = FALSE)
{
  if(is.null(names)) names <- names(height)  
  height <- as.vector(height)
  n <- length(height)
  width <- rep(width, length.out = n)
  offset <- rep(offset, length.out = n)

  if(is.null(names)) names <- rep("", n)  
  if(is.null(xlim)) xlim <- c(1 - mean(width[c(1, n)]), n + mean(width[c(1, n)]))
  if(is.null(ylim)) ylim <- c(min(offset), max(height + offset))

  if(newpage) grid.newpage()
  pushViewport(plotViewport(xscale = xlim, yscale = ylim, default.units = "native", name = name))
  grid.rect(x = 1:n, y = offset, width = width, height = height,
    just = c("centre", "bottom"), default.units = "native", gp = gp)
  grid.yaxis()
  grid.text(names, x = unit(1:n, "native"), y = unit(rep(-1.5, n), "lines"))
  grid.text(xlab, y = unit(-3.5, "lines"))
  grid.text(ylab, x = unit(-3, "lines"), rot = 90)
  grid.text(main, y = unit(1, "npc") + unit(2, "lines"), gp = gpar(fontface = "bold"))
  if(pop) popViewport() else upViewport()
  invisible(1:n)
}
