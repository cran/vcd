#################################################################333
## assocplot

assoc <- function(x, ...)
  UseMethod("assoc")

assoc.formula <-
function(formula, data = NULL, subset, na.action, ..., main = NULL)
{
    if (is.logical(main) && main)
      main <- deparse(substitute(data))

    assoc.default(ftable(formula, data, subset, na.action), main = main, ...)
}

assoc.default <- function(x,
                          row_vars = NULL, col_vars = NULL,
                          compress = TRUE, xlim = NULL, ylim = NULL,
                          spacing = spacing_conditional(sp = 0),
                          spacing_args = list(),
                          split_vertical = NULL,
                          keep_aspect_ratio = FALSE,
			  residuals_type = "Pearson",
                          xscale = 0.9, yspace = unit(0.5, "lines"),
                          main = NULL,
                          ...,
                          gp_axis = gpar(lty = 3)
                          ) {

  if (is.logical(main) && main)
    main <- deparse(substitute(x))

  if (!inherits(x, "ftable"))
    x <- structable(x)

#     {
#     if (is.null(row_vars) && is.null(col_vars) && is.table(x))
#       row_vars <- names(dimnames(x))[seq(1, length(dim(x)), by = 2)]
#     x <- ftable(x, row.vars = row_vars, col.vars = col_vars)
#   }

  tab <- as.table(x)
  dl <- length(dim(tab))
  
  ## spacing
  cond <- rep(TRUE, dl)
  cond[length(attr(x, "row.vars")) + c(0, length(attr(x, "col.vars")))] <- FALSE
  if (inherits(spacing, "panel_generator"))
    spacing <- do.call("spacing", spacing_args)
  spacing <- spacing(dim(tab), condvars = which(cond))

  ## splitting arguments
  if (is.null(split_vertical))
    split_vertical <- attr(x, "split_vertical")
  if (is.null(split_vertical)) {
    split_vertical <- rep(FALSE, dl)
    names(split_vertical) <- names(dimnames(tab))
    split_vertical[names(attr(x, "col.vars"))] <- TRUE
  }

  if(match.arg(tolower(residuals_type), "pearson") != "pearson")
    warning("Only Pearson residuals can be visualized with association plots.")
  
  strucplot(tab,
            spacing = spacing,
            split_vertical = split_vertical,
            panel = struc_assoc(compress = compress, xlim = xlim, ylim = ylim,
              yspace = yspace, xscale = xscale, gp_axis = gp_axis),
            keep_aspect_ratio = keep_aspect_ratio,
            residuals_type = "Pearson",
            main = main, 
            ...)
}

struc_assoc <- function(compress = TRUE, xlim = NULL, ylim = NULL,
                        yspace = unit(0.5, "lines"), xscale = 0.9,
                        gp_axis = gpar(lty = 3))
  function(residuals, observed = NULL, expected, spacing, gp, split_vertical) {
    dn <- dimnames(expected)
    dnn <- names(dn)
    dx <- dim(expected)
    dl <- length(dx)

    ## axis limits
    resid <- structable(residuals, split_vertical = split_vertical)
    sexpected <- structable(sqrt(expected), split_vertical = split_vertical)
    rfunc <- function(x) c(min(x, 0), max(x, 0))
    if (is.null(ylim))
      ylim <- if (compress)
        matrix(apply(resid, 1, rfunc), nrow = 2)
      else
        rfunc(resid)
    if (!is.matrix(ylim))
      ylim <- matrix(ylim, nrow = 2, ncol = nrow(resid))

    attr(ylim, "split_vertical") <- rep(TRUE, sum(!split_vertical))
    attr(ylim, "dnames") <- dn[!split_vertical]
    class(ylim) <- "structable"

    if(is.null(xlim))
      xlim <- if (compress)
        matrix(c(-0.5, 0.5) %o% apply(sexpected, 2, max), nrow = 2)
      else
        c(-0.5, 0.5) * max(sexpected)
    if (!is.matrix(xlim))
      xlim <- matrix(xlim, nrow = 2, ncol = ncol(resid))
    attr(xlim, "split_vertical") <- rep(TRUE, sum(split_vertical))
    attr(xlim, "dnames") <- dn[split_vertical]
    class(xlim) <- "structable"

    ## split workhorse
    split <- function(res, sexp, i, name, row, col) {
      v <- split_vertical[i]
      splitbase <- if (v) sexp else res
      splittab <- lapply(seq(dx[i]), function(j) splitbase[[j]])
      len <- sapply(splittab, function(x) sum(x[1,] - x[2,]))

      d <- dx[i]

      ## compute total cols/rows and build split layout
      dist <- unit.c(unit(len, "null"), spacing[[i]]  + (1 * !v) * yspace)
      idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
      layout <- if (v)
        grid.layout(ncol = 2 * d - 1, widths = dist[idx])
      else
        grid.layout(nrow = 2 * d - 1, heights = dist[idx])
      vproot <- viewport(layout.pos.col = col, layout.pos.row = row,
                         layout = layout, name = substr(name, 1, nchar(name) - 1))

      ## next level: either create further splits, or final viewports
      name <- paste(name, dnn[i], "=", dn[[i]], ",", sep = "")
      rows <- cols <- rep.int(1, d)
      if (v) cols <- 2 * 1:d - 1 else rows <- 2 * 1:d - 1

      f <- if (i < dl) {
        if (v)
          function(m) split(res, splittab[[m]], i + 1, name[m], rows[m], cols[m])
        else
          function(m) split(splittab[[m]], sexp, i + 1, name[m], rows[m], cols[m])
      } else {
        if (v)
          function(m) viewport(layout.pos.col = cols[m], layout.pos.row = rows[m],
                               name = substr(name[m], 1, nchar(name[m]) - 1),
                               yscale = res[,1],
                               xscale = sexp[,m], default.units = "null")
        else
          function(m) viewport(layout.pos.col = cols[m], layout.pos.row = rows[m],
                               name = substr(name[m], 1, nchar(name[m]) - 1),
                               yscale = res[,m],
                               xscale = sexp[,1], default.units = "null")
      }
      vpleaves <- structure(lapply(1:d, f), class = c("vpList", "viewport"))

      vpTree(vproot, vpleaves)
    }

    ## start spltting on top, creates viewport-tree
    pushViewport(split(ylim, xlim, i = 1, name = "cell:", row = 1, col = 1))

    ## draw tiles
    mnames <- paste(apply(expand.grid(dn), 1,
                          function(i) paste(dnn, i, collapse = ",", sep = "=")
                          )
                    )
    for (i in seq(along = mnames)) {
      seekViewport(paste("cell:", mnames[i], sep = ""))
      grid.lines(y = unit(0, "native"), gp = gp_axis)
      grid.rect(y = 0, x = 0,
                height = residuals[i],
                width = xscale * unit(sqrt(expected[i]), "native"),
                default.units = "native",
                gp = structure(lapply(gp, function(x) x[i]), class = "gpar"),
                just = c("center", "bottom"),
                name = paste("rect:", mnames[i], sep = "")
                )
    }

  }
class(struc_assoc) <- "panel_generator"
