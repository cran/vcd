#########################################
## structable

structable <- function(x, ...)
  UseMethod("structable")

structable.formula <- function(formula, data = NULL, direction = NULL,
                               split_vertical = FALSE, ..., subset, na.action) {
    if (missing(formula) || !inherits(formula, "formula")) 
        stop("formula is incorrect or missing")

    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    if (!is.null(direction))
      split_vertical <- direction == "v"
        
    ## only rhs present without `.' in lhs => xtabs-interface
    if (length(formula) != 3) {
      if (formula[[1]] == "~") {
        if (inherits(edata, "ftable") || inherits(edata, "table") || 
            length(dim(edata)) > 2) {
          data <- as.table(data)
          varnames <- attr(terms(formula), "term.labels")
          dnames <- names(dimnames(data))
          di <- match(varnames, dnames)
          if (any(is.na(di))) 
            stop("incorrect variable names in formula")
          if (all(varnames != ".")) 
            data <- margin.table(data, di)
          return(structable(data, split_vertical = split_vertical, ...))
        }
        else if (is.data.frame(data)) {
          if ("Freq" %in% colnames(data))
            return(structable(xtabs(formula(paste("Freq", deparse(formula))),
                                    data),
                              split_vertical = split_vertical, ...))
          else
            return(structable(xtabs(formula, data),  split_vertical = split_vertical, ...))
            
        } else {
          if (is.matrix(edata)) 
            m$data <- as.data.frame(data)
          m$... <- NULL
          m[[1]] <- as.name("model.frame")
          mf <- eval(m, parent.frame())
          return(structable(table(mf), split_vertical = split_vertical, ...))
        }

     } else
        stop("formula must have both left and right hand sides")
    }

    ## `ftable' behavior
    if (any(attr(terms(formula), "order") > 1)) 
        stop("interactions are not allowed")
    rvars <- attr(terms(formula[-2]), "term.labels")
    cvars <- attr(terms(formula[-3]), "term.labels")
    rhs.has.dot <- any(rvars == ".")
    lhs.has.dot <- any(cvars == ".")
    if (lhs.has.dot && rhs.has.dot) 
        stop(paste("formula has", sQuote("."), "in both left and right hand side"))
    if (inherits(edata, "ftable") || inherits(edata, "table") || 
        length(dim(edata)) > 2) {
        if (inherits(edata, "ftable"))
            data <- as.table(data)
        
        dnames <- names(dimnames(data))
        rvars <- pmatch(rvars, dnames)
        cvars <- pmatch(cvars, dnames)
        if (rhs.has.dot) 
            rvars <- seq(along = dnames)[-cvars]
        else if (any(is.na(rvars))) 
          stop("incorrect variable names in rhs of formula")
        if (lhs.has.dot) 
            cvars <- seq(along = dnames)[-rvars]
        else if (any(is.na(cvars))) 
          stop("incorrect variable names in lhs of formula")
        split_vertical <- c(rep(FALSE, length(rvars)), rep(TRUE, length(cvars)))
        structable(margin.table(data, c(rvars, cvars)), split_vertical = split_vertical, ...)
    } else {
        if (is.matrix(edata)) 
            m$data <- as.data.frame(data)
        m$... <- NULL
        if (!is.null(data) && is.environment(data)) {
            dnames <- names(data)
            if (rhs.has.dot) 
                rvars <- seq(along = dnames)[-cvars]
            if (lhs.has.dot) 
                cvars <- seq(along = dnames)[-rvars]
        }
        else {
            if (lhs.has.dot || rhs.has.dot) 
                stop("cannot use dots in formula with given data")
        }
        if ("Freq" %in% colnames(m$data))
          m$formula <- formula(paste("Freq~", paste(c(rvars, cvars), collapse = "+")))
        else
          m$formula <- formula(paste("~", paste(c(rvars, cvars), collapse = "+")))
        m[[1]] <- as.name("xtabs")
        mf <- eval(m, parent.frame())
        split_vertical <- c(rep(FALSE, length(rvars)), rep(TRUE, length(cvars)))
        structable(mf, split_vertical = split_vertical, ...)
    }
}
  

structable.default <- function(..., direction = NULL, split_vertical = FALSE) {
  ## several checks & transformations for arguments
  args <- list(...)
  
  if (length(args) == 0) 
    stop("Nothing to tabulate")
  
  x <- args[[1]]
  x <- if (is.list(x)) 
    table(x)
  else if (inherits(x, "ftable")) 
    as.table(x)
  else if (!(is.array(x) && length(dim(x)) > 1 || inherits(x, "table")))
    do.call("table", as.list(substitute(list(...)))[-1])
  else
    x

  ## splitting argument
  dl <- length(dim(x))
  if (!is.null(direction))
    split_vertical <- direction == "v"
  if (length(split_vertical) == 1)
    split_vertical <- rep(c(split_vertical, !split_vertical), length.out = dl)
  if (length(split_vertical) < dl)
    split_vertical <- rep(split_vertical, length.out = dl)
    

  ## permute & reshape
  ret <- aperm(x, c(rev(which(!split_vertical)), rev(which(split_vertical))))

  dn <- dimnames(x)
  rv <- dn[split_vertical]
  cv <- dn[!split_vertical]
  rl <- if (length(rv)) sapply(rv, length) else 1
  cl <- if (length(cv)) sapply(cv, length) else 1
  dim(ret) <- c(prod(cl), prod(rl))
  
  ## add dimnames
  attr(ret, "dnames") <- dn
  attr(ret, "split_vertical") <- split_vertical

  ## add dimension attributes in ftable-format
  attr(ret, "col.vars") <- rv
  attr(ret, "row.vars") <- cv

  class(ret) <- c("structable", "ftable")
  ret
}

"[[.structable" <- function(x, ...) {
  ## handle one-arg cases
  if (nargs() < 3)
    if (length(..1) > 1)
      ## resolve calls like x[[c(1,2)]]
      return(x[[ ..1[1] ]] [[ ..1[-1] ]])
    else
      ## resolve x[[foo]] 
      return(if (attr(x, "split_vertical")[1]) x[[,..1]] else x[[..1,]])

  ## handle calls like x[[c(1,2), c(3,4)]]
  if (length(..1) > 1 && length(..2) > 1)
    return(x[[ ..1[1], ..2[[1]] ]] [[ ..1[-1], ..2[-1] ]])
  ## handle calls like x[[c(1,2), 3]]
  if (length(..1) > 1)
    return(x[[ ..1[1], ..2 ]] [[ ..1[-1], ]])
  ## handle calls like x[[1, c(1,3)]]
  if (length(..2) > 1)
    return(x[[ ..1, ..2[[1]] ]] [[ , ..2[-1] ]])

  ## final cases like x[[1,2]] or x[[1,]] or x[[,1]]
  rv <- attr(x, "dnames")[!attr(x, "split_vertical")]
  cv <- attr(x, "dnames")[attr(x, "split_vertical")]

  if (!is.symbol(..1)) rstep <- dim(x)[1] / length(rv[[1]])
  if (!is.symbol(..2)) cstep <- dim(x)[2] / length(cv[[1]])

  ret <- x[if (!is.symbol(..1)) (1 + (..1 - 1) * rstep) : (..1 * rstep) else 1:nrow(x),
           if (!is.symbol(..2)) (1 + (..2 - 1) * cstep) : (..2 * cstep) else 1:ncol(x),
           drop = FALSE
           ]

  attr(ret, "split_vertical") <- attr(x, "split_vertical")
  attr(ret, "dnames") <- attr(x, "dnames")
  
  if (!is.symbol(..1)) {
    i <- which(!attr(ret, "split_vertical"))[1]
    attr(ret, "split_vertical") <- attr(ret, "split_vertical")[-i]
    attr(ret, "dnames") <- attr(ret, "dnames")[-i]
  }
    

  if (!is.symbol(..2)) {
    i <- which(attr(ret, "split_vertical"))[1]
    attr(ret, "split_vertical") <- attr(ret, "split_vertical")[-i]
    attr(ret, "dnames") <- attr(ret, "dnames")[-i]
  }

  class(ret) <- class(x)
  ret
}
