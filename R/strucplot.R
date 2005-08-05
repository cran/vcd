################################################################
### strucplot - generic plot framework for mosaic-like layouts
### 2 panel functions are provided: struc_mosaic and struc_assoc
################################################################

strucplot <- function(## main parameters
                      x,
                      residuals = NULL,
                      expected = NULL,
		      condvars = NULL,
                      shade = NULL,
                      type = c("observed", "expected"),
                      residuals_type = c("Pearson", "deviance", "FT"),
                      
                      ## layout
                      split_vertical = TRUE, 
                      spacing = spacing_equal,
                      spacing_args = list(),
                      gp = NULL,
		      gp_args = list(),   
                      labeling = labeling_border,
                      labeling_args = list(),
                      panel = struc_mosaic,
                      panel_args = list(),
                      legend = NULL,
                      legend_args = list(),
                      
                      main = NULL,
                      sub = NULL,
                      margins = rep.int(2.5, 4),
                      legend_width = unit(5, "lines"),
                      
                      ## control parameters
                      title_gp = gpar(fontsize = 20),
                      newpage = TRUE,
                      pop = TRUE,
                      keep_aspect_ratio = TRUE
                      ) {
  ## default behaviour of shade
  if (is.null(shade)) shade <- !is.null(gp) || !is.null(expected)
		      
  type <- match.arg(type)
  residuals_type <- match.arg(tolower(residuals_type), c("pearson", "deviance", "ft"))

  ## table characteristics
  d <- dim(x)
  dl <- length(d)
  dn <- dimnames(x)
  if (is.null(dn))
    dn <- dimnames(x) <- lapply(d, seq)
  dnn <- names(dimnames(x))
  if (is.null(dnn))
    dnn <- names(dn) <- names(dimnames(x)) <- LETTERS[1:dl]

  ## replace NAs by 0
  if (any(nas <- is.na(x))) x[nas] <- 0

  ## model fitting
  ## A parameter df is added for inference
  ## (which is done in the shading (generating) functions).
  df <- NULL
  if (is.null(expected) || !is.numeric(expected))
    if (inherits(expected, "formula")) {
      fm <- loglm(expected, x, fitted = TRUE)
      expected <- fitted(fm)
      df <- fm$df
    } else {
      if (is.null(expected))
        expected <- if (is.null(condvars))
          as.list(1:dl)
        else
          lapply((condvars + 1):dl, c, seq(condvars))

      fm <- loglin(x, expected, fit = TRUE, print = FALSE)
      expected <- fm$fit
      df <- fm$df
    }
  
  ## compute residuals
  if (is.null(residuals))
    residuals <- switch(residuals_type,
                        pearson = (x - expected) / sqrt(ifelse(expected > 0, expected, 1)),
                        deviance = {
                          tmp <- 2 * (x * log(ifelse(x == 0, 1, x / ifelse(expected > 0, expected, 1))) - (x - expected))
                          tmp <- sqrt(pmax(tmp, 0))
                          ifelse(x > expected, tmp, -tmp)
                        },
                        ft = sqrt(x) + sqrt(x + 1) - sqrt(4 * expected + 1)
                        )
  ## replace NAs by 0
  if (any(nas <- is.na(residuals))) residuals[nas] <- 0

  ## splitting
  if (length(split_vertical) == 1)
    split_vertical <- rep(c(split_vertical, !split_vertical), length.out = dl)

  ## spacing
  if (is.function(spacing)) {
    if (inherits(spacing, "panel_generator"))
      spacing <- do.call("spacing", spacing_args)
    spacing <- spacing(d, condvars)
  }

  ## gp (color, fill, lty, etc.) argument
  if (shade) {
    if (is.null(gp)) gp <- shading_hcl
    if (is.function(gp)) {
      if (is.null(legend) || (is.logical(legend) && legend))
        legend <- legend_resbased
      gpfun <- if(inherits(gp, "panel_generator"))
        do.call("gp", c(list(x, residuals, expected, df), as.list(gp_args))) else gp
      gp <- gpfun(residuals)
    } else if (!is.null(legend) && !(is.logical(legend) && !legend))
      stop("gp argument must be a shading function for drawing a legend")
  } else {
    if(!is.null(gp)) {
      warning("gp parameter ignored since shade=FALSE")
      gp <- NULL
    }
  }

  ## choose gray when no shading is used
  if (is.null(gp)) gp <- gpar(fill = grey(0.8))

  ## recycle gpar values in the last dimension
  size <- prod(d)
  FUN <- function(par) if (length(par) < size) aperm(array(par, dim = rev(d))) else par
  gp <- structure(lapply(gp, FUN), class = "gpar")
  
  ## set up page
  if (newpage)
    grid.newpage()
  if (keep_aspect_ratio)
    pushViewport(viewport(width = 1, height = 1, default.unit = "snpc"))
  
  pushViewport(vcdViewport(mar = margins,
                           legend = shade && !(is.null(legend) || is.logical(legend) && !legend),
                           main = !is.null(main), sub = !is.null(sub), keep_aspect_ratio = keep_aspect_ratio,
                           legend_width = legend_width))

  ## legend
  if (inherits(legend, "panel_generator"))
    legend <- do.call("legend", legend_args)
  if (shade && !is.null(legend) && !(is.logical(legend) && !legend)) {
    seekViewport("legend")
    residuals_type = switch(residuals_type, deviance = "deviance", ft = "Freeman-Tukey", pearson = "Pearson")
    legend(residuals, gpfun, paste(residuals_type, "residuals:", sep = "\n"))
  }

  ## titles
  if (!is.null(main)) {
    seekViewport("main")
    if (is.logical(main) && main)
      main <- deparse(substitute(x))
    grid.text(main, gp = title_gp)
  }

  if (!is.null(sub)) {
    seekViewport("sub")
    grid.text(sub, gp = title_gp)
  }

  ## make plot
  seekViewport("plot")
  
  if (inherits(panel, "panel_generator"))
    panel <- do.call("panel", panel_args)
  panel(residuals = residuals,
        observed = if (type == "observed") x else expected,
        expected = if (type == "observed") expected else x,
        spacing = spacing,
        gp = gp,
        split_vertical = split_vertical)

  upViewport(dl)

  ## labels
  if (!is.null(labeling)) {
    if (inherits(labeling, "panel_generator"))
      labeling <- do.call("labeling", labeling_args)
    labeling(dn, split_vertical, condvars)
  }

  ## pop/move up viewport

  seekViewport("base") 
  ## one more up if sandwich-mode
  if (!is.null(main) || !is.null(sub) ##||
##      (shade && !is.null(legend) && !(is.logical(legend) && !legend))
      ) upViewport()
  if (pop) popViewport(1 + keep_aspect_ratio) else upViewport(1 + keep_aspect_ratio)

  ## return visualized table
  invisible(structable(if (type == "observed") x else expected,
                      split_vertical = split_vertical))
}

vcdViewport <- function(mar = rep.int(2.5, 4),
                        legend_width = unit(5, "lines"),
                        legend = FALSE, main = FALSE, sub = FALSE,
                        keep_aspect_ratio = TRUE)
{
  mar <- if (!is.unit(mar))
    unit(pexpand(mar, 4, rep.int(2.5, 4), c("top","right","bottom","left")), "lines")
  else
    unit.rep(mar, length.out = 4)
  if (!is.unit(legend_width))
    legend_width <- unit(legend_width, "lines")
  vpPlot <- vpStack(viewport(layout.pos.col = 2, layout.pos.row = 2),
                    viewport(width = 1, height = 1, name = "plot",
                             default.units = if (keep_aspect_ratio) "snpc" else "npc"))
  vpMarginBottom <- viewport(layout.pos.col = 2, layout.pos.row = 3, name = "margin_bottom")
  vpMarginLeft <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "margin_left")
  vpMarginTop <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "margin_top")
  vpMarginRight <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "margin_right")
  vpCornerTL <- viewport(layout.pos.col = 1, layout.pos.row = 1, name = "corner_top_left")
  vpCornerTR <- viewport(layout.pos.col = 3, layout.pos.row = 1, name = "corner_top_right")
  vpCornerBL <- viewport(layout.pos.col = 1, layout.pos.row = 3, name = "corner_bottom_left")
  vpCornerBR <- viewport(layout.pos.col = 3, layout.pos.row = 3, name = "corner_bottom_right")

  if (legend) {
    vpLegend <- viewport(layout.pos.col = 4, layout.pos.row = 2, name = "legend")
    vpPval <- viewport(layout.pos.col = 4, layout.pos.row = 3, name = "pval")
    vpBase <- viewport(layout.pos.row = 1 + main,
                       layout = grid.layout(3, 4,
                         widths = unit.c(mar[4], unit(1, "null"), mar[2], legend_width),
                         heights = unit.c(mar[1], unit(1, "null"), mar[3])),
                       name = "base")
    vpPlotregion <- vpTree(vpBase, vpList(vpMarginBottom, vpMarginLeft, vpMarginTop,
                                          vpMarginRight, vpPval, vpLegend,
                                          vpCornerTL, vpCornerTR, vpCornerBL,
                                          vpCornerBR, vpPlot))
  } else {
    vpBase <- viewport(layout.pos.row = 1 + main,
                       layout = grid.layout(3, 3,
                         widths = unit.c(mar[4], unit(1, "null"), mar[2]),
                         heights = unit.c(mar[1], unit(1, "null"), mar[3])
                         ),
                       name = "base")
    vpPlotregion <- vpTree(vpBase,
                           vpList(vpMarginBottom, vpMarginLeft, vpMarginTop, vpMarginRight,
                                  vpCornerTL, vpCornerTR, vpCornerBL, vpCornerBR, vpPlot))
  }

  ## main/sub-title, margins for legend layout
  if (main || sub) {
    vpTop <- viewport(layout.pos.row = 1, name = "main")
    vpSub <- viewport(layout.pos.row = 2 + main, name = "sub")
    
    sandwich <-
## no additional space when keep_aspect_ratio = F
#       if (legend) {
#       space <- max(legend_width + mar[2] + mar[4] - mar[1] - mar[3],
#                    unit((main + sub) * 2, "lines"))
#       vplist <- vpList(vpTop, vpPlotregion, vpSub)
#       viewport(layout = grid.layout(3, 1,
#                  height = unit.c(0.5 * space, unit(1, "null"), 0.5 * space)))
#     } else
    if (main && sub) {
      vplist <- vpList(vpTop, vpPlotregion, vpSub)
      viewport(layout = grid.layout(3, 1,
                 height = unit.c(unit(2, "lines"), unit(1, "null"), unit(2, "lines"))))
    } else if (main) {
      vplist <- vpList(vpTop, vpPlotregion)
      viewport(layout = grid.layout(2, 1,
                 height = unit.c(unit(2, "lines"), unit(1, "null"))))
    } else {
      vplist <- vpList(vpPlotregion, vpSub)
      viewport(layout = grid.layout(2, 1,
                 height = unit.c(unit(1, "null"), unit(2, "lines"))))
    }

    vpTree(sandwich, vplist)
  } else vpPlotregion
}
