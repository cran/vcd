cotabplot <- function(x, ...)
{
  UseMethod("cotabplot")
}

cotabplot.formula <- function(formula, data = NULL, ...)
{
  m <- match.call()
  edata <- eval(m$data, parent.frame())
  
  fstr <- deparse(formula)
  fstr <- gsub("*", "+", fstr, extended = FALSE, fixed = TRUE)
  fstr <- gsub("/", "+", fstr, extended = FALSE, fixed = TRUE)
  fstr <- gsub("(", "", fstr, extended = FALSE, fixed = TRUE)
  fstr <- gsub(")", "", fstr, extended = FALSE, fixed = TRUE)  
  
  fstr <- strsplit(paste(fstr, collapse = ""), "~")
  vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
  varnames <- vars[[1]]
  condnames <- if(length(vars) > 1) vars[[2]] else NULL

  if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
    tab <- as.table(data)
    if(all(varnames != ".")) {
      ind <- match(varnames, names(dimnames(tab)))
      if (any(is.na(ind)))
        stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))
      
      if (!is.null(condnames)) {
        condind <- match(condnames, names(dimnames(tab)))
        if (any(is.na(condind)))
          stop(paste("Can't find", paste(condnames[is.na(condind)], collapse=" / "), "in", deparse(substitute(data))))
        ind <- c(condind, ind)
      }
      tab <- margin.table(tab, ind)
    }
  } else {
    tab <- if ("Freq" %in% colnames(data))
      xtabs(formula(paste("Freq~", paste(c(condnames, varnames), collapse = " + "))),
            data = data)
    else
      xtabs(formula(paste("~", paste(c(condnames, varnames), collapse = " + "))),
            data = data)
  }
  tab <- margin.table(tab, match(c(varnames, condnames), names(dimnames(tab))))
  cotabplot(tab, cond = condnames, ...)
}


cotabplot.default <- function(x, cond = NULL,
  panel = cotab_mosaic, panel_args = list(),
  margins = rep(1, 4),
  text_gp = gpar(fontsize = 12), rect_gp = gpar(fill = grey(0.9)),
  pop = TRUE, newpage = TRUE,
  ...)
{
  ## coerce to table
  x <- as.table(x)
  
  ## initialize newpage
  if(newpage) grid.newpage()

  ## process default option
  ldx <- length(dim(x))
  if(is.null(cond)) {
    indep <- if(ldx > 1) 1:2 else 1
    if(ldx > 2) cond <- 3:ldx
  } else {
    if(is.character(cond)) cond <- match(cond, names(dimnames(x)))
    cond <- as.integer(cond)
    indep <- (1:ldx)[!(1:ldx %in% cond)] 
  }

  ## sort margins      
  x <- margin.table(x, c(indep, cond))

  ## convenience variables that describe conditioning variables  
  if(is.null(cond)) {
    cond.n <- 0
    cond.num <- cond.dnam <- cond.char <- NULL
  } else {
    cond.n <- length(cond)               ## number of variables
    cond.num <- (length(indep) + 1):ldx  ## position in x
    cond.dnam <- dimnames(x)[cond.num]   ## corresponding dimnames
    cond.char <- names(cond.dnam)        ## names of variables
  }

  ## create panel function (if necessary)
  if(inherits(panel, "grapcon_generator"))
    panel <- do.call("panel", c(list(x, cond.char), as.list(panel_args), list(...)))

  if(cond.n < 1) panel(x, NULL) ## no conditioning variables
  else {
  
  cond.nlevels <- sapply(cond.dnam, length)
  nplots <- prod(cond.nlevels)
  condition <- as.matrix(expand.grid(cond.dnam))

  ## compute layout
  #Z# needs fixing for more than two conditioning variables
  layout <- c(1,1,1) ## rows, cols, pages
  if(cond.n == 1) {
    layout[2] <- ceiling(sqrt(floor(cond.nlevels)))
    layout[1] <- ceiling(cond.nlevels/layout[2])
    layout <- expand.grid(lapply(layout, function(x) 1:x))[1:nplots,]
  }
  else {
    layout[1] <- cond.nlevels[1]
    layout[2] <- cond.nlevels[2]
    if(cond.n >= 3) layout[3] <- nplots/prod(cond.nlevels[1:2]) #Z# FIXME
    if(layout[3] > 1) stop("multiple pages not supported yet")
    layout <- expand.grid(lapply(layout, function(x) 1:x))
  }

  ## push basic grid of nr x nc cells
  nr <- max(layout[,1])
  nc <- max(layout[,2])
  pushViewport(plotViewport(margins))
  pushViewport(viewport(layout = grid.layout(nr, nc, widths = unit(1/nc, "npc"))))

  strUnit <- unit(2 * ncol(condition), "strheight", "A")
  cellport <- function(name) viewport(layout = grid.layout(2, 1,
        heights = unit.c(strUnit, unit(1, "npc") - strUnit)),
	name = name)

  ## go through each conditioning combination
  for(i in 1:nrow(condition)) {

    ## conditioning information in ith cycle
    condi <- as.vector(condition[i,])
    names(condi) <- colnames(condition)
    condistr <- paste(condi, collapse = ".")
    condilab <- paste(cond.char, condi, sep = " = ")

    ## header
    pushViewport(viewport(layout.pos.row = layout[i,1], layout.pos.col = layout[i,2]))
    pushViewport(cellport(paste("cell", condistr, sep = ".")))
    pushViewport(viewport(layout.pos.row = 1, name = paste("lab", condistr, sep = ".")))
    grid.rect(gp = rect_gp)
    grid.text(condilab, y = cond.n:1/cond.n - 1/(2*cond.n), gp = text_gp)
    grid.segments(0, 0:cond.n/cond.n, 1, 0:cond.n/cond.n)
    upViewport()

    ## main plot
    pushViewport(viewport(layout.pos.row = 2, name = paste("plot", condistr, sep = ".")))
    panel(x, condi)
    upViewport(2)
    grid.rect(gp = gpar(fill = "transparent"))
    upViewport()
  }
  upViewport()
  if(pop) popViewport() else upViewport()
  }

  invisible(x)
}

cotab_mosaic <- function(x = NULL, condvars = NULL, ...) {
  function(x, condlevels) {
    if(is.null(condlevels)) mosaic(x, newpage = FALSE, pop = TRUE, ...)
      else mosaic(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
                  newpage = FALSE, pop = TRUE, ...)
  }
}
class(cotab_mosaic) <- "grapcon_generator"

cotab_sieve <- function(x = NULL, condvars = NULL, ...) {
  function(x, condlevels) {
    if(is.null(condlevels)) sieve(x, newpage = FALSE, pop = TRUE, ...)
      else sieve(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
                 newpage = FALSE, pop = TRUE, ...)
  }
}
class(cotab_sieve) <- "grapcon_generator"

cotab_assoc <- function(x = NULL, condvars = NULL, ylim = NULL, ...) {
  if(!is.null(x)) {
    fm <- coindep_test(x, condvars, n = 1)
    if(is.null(ylim)) ylim <- range(residuals(fm))
  }
  
  function(x, condlevels) {
    if(is.null(condlevels)) assoc(x, newpage = FALSE, pop = TRUE, ylim = ylim, ...)
      else assoc(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
                  newpage = FALSE, pop = TRUE, ylim = ylim, ...)
  }
}
class(cotab_assoc) <- "grapcon_generator"

cotab_fourfold <- function (x = NULL, condvars = NULL, ...) {
  function(x, condlevels) {
    if (is.null(condlevels)) 
      fourfold(x, newpage = FALSE, ...)
    else
      fourfold(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
               newpage = FALSE, ...)
  }
}
class(cotab_fourfold) <- "grapcon_generator"

cotab_coindep <- function(x, condvars,
  test = c("max", "Chisq"), level = NULL, n = 1000,
  h = NULL, c = NULL, l = NULL, lty = 1,
  type = c("mosaic", "assoc"),
  legend = FALSE, ylim = NULL, ...)
{
  if(is.null(condvars))
    stop("at least one conditioning variable is required")

  ## set color defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- c(100, 20)
  if(is.null(l)) l <- c(90, 50)  
  
  ## process conditional variables and get independent variables
  ## store some convenience information
  ldx <- length(dim(x))
  if(is.character(condvars)) condvars <- match(condvars, names(dimnames(x)))
  condvars <- as.integer(condvars)
  indep <- (1:ldx)[!(1:ldx %in% condvars)] 
  
  ## sort margins      
  x <- margin.table(x, c(indep, condvars))

  ind.n <- length(indep)            
  ind.num <- 1:ind.n                
  ind.dnam <- dimnames(x)[ind.num]  
  ind.char <- names(ind.dnam)       
  cond.n <- length(condvars)		
  cond.num <- (ind.n + 1):length(dim(x))
  cond.dnam <- dimnames(x)[cond.num]	
  cond.char <- names(cond.dnam) 	

  test <- match.arg(test)
  if(test == "max") {
    if(is.null(level)) level <- c(0.9, 0.99)
    fm <- coindep_test(x, cond.num, n = n)
    resids <- residuals(fm)

    col.bins <- fm$qdist(sort(level))
    gpfun <- shading_hcl(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
      h = h, c = c, l = l, interpolate = col.bins, lty = lty, p.value = fm$p.value)
  } else {
    if(is.null(level)) level <- 0.95
    level <- level[1]
    
    indep.formula <- as.formula(paste("~ (",
                                      paste(ind.char, collapse = " + "),
	    			      ") * ",
				      paste(cond.char, collapse = " * "),
				      sep = ""))
    fm <- loglm(indep.formula, data = x, fitted = TRUE)
    resids <- residuals(fm, type = "pearson")
    
    gpfun <- shading_hcl(x, fitted(fm), resids, fm$df,
      h = h, c = c, l = l, lty = lty)
  }

  ## if "maxChisq" should be implemented, use something
  ## along the lines of
  ## mosaiclist <- list()
  ## 
  ## rval <- function(x, indep, cond)
  ##   mosaiclist[[ ... ]](co_table(margin.table(x, c(indep, cond)), cond)[[ ... ]])

  type <- match.arg(type)
  if(type == "mosaic") {
    rval <- function(x, condlevels) {
    if(is.null(condlevels))
      mosaic(x, newpage = FALSE, pop = TRUE,
             gp = gpfun, legend = legend, ...)
    else
      mosaic(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
             newpage = FALSE, pop = TRUE,
	     gp = gpfun, legend = legend, ...)
    }
  } else {
    if(is.null(ylim)) ylim <- range(resids)

    rval <- function(x, condlevels) {
    if(is.null(condlevels))
      assoc(x, newpage = FALSE, pop = TRUE,
             gp = gpfun, legend = legend, ylim = ylim, ...)
    else
      assoc(co_table(x, names(condlevels))[[paste(condlevels, collapse = ".")]],
             newpage = FALSE, pop = TRUE,
	     gp = gpfun, legend = legend, ylim = ylim, ...)
    }
  }

  return(rval)
}
class(cotab_coindep) <- "grapcon_generator"


## Examples:
## these work (although they are not very meaningful)
## 
## cotabplot(~ Class + Survived | Age + Sex, data = Titanic)
## cotabplot(Titanic)
## cotabplot(Titanic, panel = cotab_assoc)

## here are some bugs in legends function
## 
## cotabplot(Titanic, shade = TRUE)
## cotabplot(Titanic, shade = TRUE, legend = FALSE)

## further illustrations (legends need to be switched on, see above)
## 
## cotabplot(~ Gender + Admit | Dept, data = UCBAdmissions)
## cotabplot(~ Gender + Admit | Dept, data = UCBAdmissions,
##   panel = cotab_coindep, legend = TRUE)
## cotabplot(~ Gender + Admit | Dept, data = UCBAdmissions,
##   panel = cotab_coindep, legend = TRUE, type = "assoc")
