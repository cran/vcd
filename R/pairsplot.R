#################################################################
### pairsplot

pairs.table <- function(x,
                  upper_panel = pairs_mosaic,
                  upper_panel_args = list(),
                        
                  lower_panel = pairs_mosaic,
                  lower_panel_args = list(),
                        
                  diag_panel = pairs_barplot,
                  diag_panel_args = list(),
                  
                  main = NULL,
                  title_gp = gpar(fontsize = 20),
                  
                  space = 0.1,
                  newpage = TRUE,
                  ...)
{
  if (newpage) grid.newpage()

  if (inherits(upper_panel, "panel_generator"))
    upper_panel <- do.call("upper_panel", upper_panel_args)
  if (inherits(lower_panel, "panel_generator"))
    lower_panel <- do.call("lower_panel", lower_panel_args)
  if (inherits(diag_panel, "panel_generator"))
    diag_panel <- do.call("diag_panel", diag_panel_args)
  
  d <- length(dim(x))
  l <- grid.layout(d, d)
  pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc")))
  pushViewport(viewport(layout = l,
                         height = if (is.null(main)) 1 else 0.9,
                         y = 0, just = "bottom"))

  if (is.logical(main) && main)
    main <- deparse(substitute(x))
  grid.text(main, y = unit(1.05, "npc"), gp = title_gp)


  for (i in 1:d)
    for(j in 1:d) {
      pushViewport(viewport(layout.pos.col = i, layout.pos.row = j))
      pushViewport(viewport(width = 1 - space, height = 1 - space))

      if (i > j)
        upper_panel(x, j, i)
      else if (i < j)
        lower_panel(x, j, i)
      else
        diag_panel(x, i)

      popViewport(2)
    }
  popViewport(2)
  invisible(x)
}

## upper/lower panels

pairs_assoc <- function(...) pairs_strucplot(panel = assoc, ...)
class(pairs_assoc) <- "panel_generator"

pairs_mosaic <- function(...) pairs_strucplot(panel = mosaic, ...)
class(pairs_mosaic) <- "panel_generator"

pairs_strucplot <- function(panel = mosaic,
                            type = c("pairwise", "total", "conditional", "joint"),
                            legend = FALSE, margins = c(0, 0, 0, 0),
                            labeling = NULL, shade = FALSE, ...) {
  type = match.arg(type)
  function(x, i, j) {
    index <- 1:length(dim(x))
    rest <- index[!index %in% c(i, j)]
    panel(x = margin.table(x, if (type == "pairwise") c(j, i) else c(j, i, rest)),
           expected = switch(type,
             pairwise =, total = NULL,
             conditional = list(c(j, rest), c(i, rest)),
             joint = list(c(j, i), rest)
             ),
           condvars = if (type == "conditional") rest else NULL,
           
           labeling = labeling,
           margins = margins,
           legend = legend,
           shade = shade,

           split_vertical = TRUE,
           
           newpage = FALSE,
           pop = TRUE,
           ...)
  }
}
class(pairs_strucplot) <- "panel_generator"
  
## diagonal panels

pairs_text <- function(dimnames = TRUE,
                       gp_vartext = gpar(fontsize = 17),
                       gp_leveltext = gpar(),
                       gp_border = gpar(),
                       ...)

  function(x, i) {
    x <- margin.table(x, i)
    grid.rect(gp = gp_border)
    grid.text(names(dimnames(x)), gp = gp_vartext,  y = 0.5 + dimnames * 0.05, ...)
    if (dimnames)
      grid.text(paste("(",paste(names(x), collapse = ","), ")", sep = ""),
                y = 0.4, gp = gp_leveltext)
  }
class(pairs_text) <- "panel_generator"

pairs_barplot <- function(gp_bars = gpar(fill = "gray"),
                          gp_vartext = gpar(fontsize = 17),
                          gp_leveltext = gpar(),
                          ...)
  function(x, i) {
    x <- margin.table(x, i)
    pushViewport(viewport(x = 0.3, y = 0.1, width = 0.7, height = 0.7,
                          yscale = c(0,max(x)), just = c("left", "bottom"))
                 )
    xpos <- seq(0, 1, length = length(x) + 1)[-1]
    halfstep <- (xpos[2] - xpos[1]) / 2
    grid.rect(xpos - halfstep, rep.int(0, length(x)), height = x,
              just = c("center", "bottom"), width = halfstep,
              gp = gp_bars, default = "native", ...)
    grid.yaxis(at = pretty(c(0,max(x))))
    grid.text(names(x), y = unit(-0.15, "npc"),
              x = xpos - halfstep, just = c("center", "bottom"), gp = gp_leveltext)
    popViewport(1)
    grid.text(names(dimnames(x)), y = 1, just = c("center", "top"), gp = gp_vartext)

  }
class(pairs_barplot) <- "panel_generator"


