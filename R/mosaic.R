###########################################################
## mosaicplot

mosaic <- function(x, ...)
  UseMethod("mosaic")

mosaic.formula <-
function(formula, data = NULL, ..., main = NULL, sub = NULL, subset = NULL)
{
  if (is.logical(main) && main)
    main <- deparse(substitute(data))
  else if (is.logical(sub) && sub)
    sub <- deparse(substitute(data))

  m <- match.call(expand.dots = FALSE)
  edata <- eval(m$data, parent.frame())
  
  fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
  vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
  varnames <- vars[[1]]
  condnames <- if (length(vars) > 1) vars[[2]] else NULL

  if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
    condind <- NULL
    dat <- as.table(data)
    if(all(varnames != ".")) {
      ind <- match(varnames, names(dimnames(dat)))
      if (any(is.na(ind)))
        stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))
      
      if (!is.null(condnames)) {
        condind <- match(condnames, names(dimnames(dat)))
        if (any(is.na(condind)))
          stop(paste("Can't find", paste(condnames[is.na(condind)], collapse=" / "), "in", deparse(substitute(data))))
        ind <- c(condind, ind)
      }
      dat <- margin.table(dat, ind)
    }
    mosaic.default(dat, main = main, sub = sub,
                   condvars = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
  } else {
      m <- m[c(1, match(c("formula", "data", "subset"), names(m), 0))]
      m[[1]] <- as.name("xtabs")
      m$formula <-
          formula(paste(if("Freq" %in% colnames(data)) "Freq",
                        "~",
                        paste(c(condnames, varnames), collapse = "+")))
      tab <- eval(m, parent.frame())
      mosaic.default(tab, main = main, sub = sub, ...)  
  }
}

mosaic.default <- function(x, condvars = NULL,
                           split_vertical = NULL, direction = NULL,
                           spacing = NULL, spacing_args = list(),
                           zero_size = 0.5, main = NULL, sub = NULL, ...) {
  if (is.logical(main) && main)
    main <- deparse(substitute(x))
  else if (is.logical(sub) && sub)
    sub <- deparse(substitute(x))

  if (is.structable(x)) {
    if (is.null(direction) && is.null(split_vertical))
      split_vertical <- attr(x, "split_vertical")
    x <- as.table(x)
  }
  if (is.null(split_vertical))
    split_vertical <- FALSE
  
  dl <- length(dim(x))

  ## splitting argument
  if (!is.null(direction))
    split_vertical <- direction == "v"
  if (length(split_vertical) == 1)
    split_vertical <- rep(c(split_vertical, !split_vertical), length.out = dl)
  if (length(split_vertical) < dl)
    split_vertical <- rep(split_vertical, length.out = dl)

  ## condvars
  if (!is.null(condvars)) {
    if (is.character(condvars))
      condvars <- match(condvars, names(dimnames(x)))
    x <- aperm(x, c(condvars, seq(dl)[-condvars]))
    if (is.null(spacing))
      spacing <- spacing_conditional
  }
  
  ## spacing argument
  if (is.null(spacing))
    spacing <- if (dl < 3) spacing_equal else spacing_increase

  strucplot(x,
            condvars = if (is.null(condvars)) NULL else length(condvars),
            core = struc_mosaic(zero_size = zero_size),
            split_vertical = split_vertical,
            spacing = spacing,
            spacing_args = spacing_args,
            main = main,
            sub = sub,
            ...)
}

struc_mosaic <- function(zero_size = 0.5)
  function(residuals, observed, expected = NULL, spacing, gp, split_vertical) {
    dn <- dimnames(observed)
    dnn <- names(dn)
    dx <- dim(observed)
    dl <- length(dx)

    ## split workhorse
    split <- function(x, i, name, row, col) {
      cotab <- co_table(x, 1)
      margin <- sapply(cotab, sum)
      v <- split_vertical[i]
      d <- dx[i]

      ## compute total cols/rows and build split layout
      dist <- unit.c(unit(margin, "null"), spacing[[i]])
      idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
      layout <- if (v)
        grid.layout(ncol = 2 * d - 1, widths = dist[idx])
      else
        grid.layout(nrow = 2 * d - 1, heights = dist[idx])
      vproot <- viewport(layout.pos.col = col, layout.pos.row = row,
                         layout = layout, name = remove_trailing_comma(name))
      
      ## next level: either create further splits, or final viewports
      name <- paste(name, dnn[i], "=", dn[[i]], ",", sep = "")
      row <- col <- rep.int(1, d)
      if (v) col <- 2 * 1:d - 1 else row <- 2 * 1:d - 1
      f <- if (i < dl) 
        function(m) split(cotab[[m]], i + 1, name[m], row[m], col[m])
      else
        function(m) viewport(layout.pos.col = col[m], layout.pos.row = row[m],
                             name = remove_trailing_comma(name[m]))
      vpleaves <- structure(lapply(1:d, f), class = c("vpList", "viewport"))

      vpTree(vproot, vpleaves)
    }

    ## start spltting on top, creates viewport-tree
    pushViewport(split(observed + .Machine$double.eps,
                       i = 1, name = "cell:", row = 1, col = 1))

    ## draw rectangles
    mnames <-  apply(expand.grid(dn), 1,
                     function(i) paste(dnn, i, collapse=",", sep = "=")
                     )
    zeros <- observed <= .Machine$double.eps

    for (i in seq(along = mnames)) {
      seekViewport(paste("cell:", mnames[i], sep = ""))
      gpobj <- structure(lapply(gp, function(x) x[i]), class = "gpar")
      if (!zeros[i]) {
        grid.rect(gp = gpobj, name = paste("rect:", mnames[i], sep = ""))
      } else { 
        grid.lines(x = 0.5, gp = gpobj)
        grid.lines(y = 0.5, gp = gpobj)
        if (zero_size > 0) {
          grid.points(0.5, 0.5, pch = 19, size = unit(zero_size, "char"),
                      gp = gpar(col = gp$fill[i]),
                      name = paste("disc:", mnames[i], sep = ""))
          grid.points(0.5, 0.5, pch = 1, size = unit(zero_size, "char"),
                      name = paste("circle:", mnames[i], sep = ""))
        }
      }
    }

  }
class(struc_mosaic) <- "grapcon_generator"
