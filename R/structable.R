#########################################
## structable

structable <- function(x, ...)
  UseMethod("structable")

structable.formula <- function(formula, data = NULL, direction = NULL,
                               split_vertical = NULL, ..., subset, na.action) {
    if (missing(formula) || !inherits(formula, "formula")) 
        stop("formula is incorrect or missing")

    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    if (!is.null(direction))
      split_vertical <- direction == "v"

    if (is.structable(data)) {
      split_vertical <- attr(data, "split_vertical")
      data <- as.table(data)
    }
    
    if (is.null(split_vertical))
      split_vertical <- FALSE
    
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
          m$... <- m$split_vertical <- m$direction <- NULL
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
        m$... <- m$split_vertical <- m$direction <- NULL
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

"[.structable" <- function(x, ...) {
  args <- if (nargs() < 3)
    list(..1)
  else
    list(..1, ..2)
  
  ## handle one-arg cases
  if (nargs() < 3)
    if (length(args[[1]]) > 1)
      ## resolve calls like x[c(1,2)]
      return(x[ args[[1]][1] ] [ args[[1]][-1] ])
    else
      ## resolve x[foo] 
      return(if (attr(x, "split_vertical")[1]) x[,args[[1]] ] else x[args[[1]],])

  ## handle calls like x[c(1,2), c(3,4)]
  if (length(args[[1]]) > 1 && length(args[[2]]) > 1)
    return(x[ args[[1]][1], args[[2]][1] ] [ args[[1]][-1], args[[2]][-1] ])
  
  ## handle calls like x[c(1,2), 3]
  if (length(args[[1]]) > 1)
    return(x[ args[[1]][1], args[[2]] ] [ args[[1]][-1], ])
  
  ## handle calls like x[1, c(1,3)]
  if (length(args[[2]]) > 1)
    return(x[ args[[1]], args[[2]][1] ] [ , args[[2]][-1] ])

  ## final cases like x[1,2] or x[1,] or x[,1]
  dnames <- attr(x, "dnames")
  split <- attr(x, "split_vertical")
  rv <- dnames[!split]
  cv <- dnames[split]

  lsym <- is.symbol(args[[1]])
  rsym <- is.symbol(args[[2]])
  if (!lsym) {
    rstep <- dim(unclass(x))[1] / length(rv[[1]])
    if (is.character(args[[1]]))
      args[[1]] <- match(args[[1]], rv[[1]])
  }
  if (!rsym) {
    cstep <- dim(unclass(x))[2] / length(cv[[1]])
    if (is.character(args[[2]]))
      args[[2]] <- match(args[[2]], cv[[1]])
  }

  lind <- if (!lsym)
    (1 + (args[[1]] - 1) * rstep) : (args[[1]] * rstep)
  else
    1:nrow(unclass(x))
  rind <- if (!rsym)
    (1 + (args[[2]] - 1) * cstep) : (args[[2]] * cstep)
  else
    1:ncol(unclass(x))
  ret <- unclass(x)[lind, rind, drop = FALSE]

  if (!lsym) {
    i <- which(!split)[1]
    split <- split[-i]
    dnames <- dnames[-i]
  }
    
  if (!rsym) {
    i <- which(split)[1]
    split <- split[-i]
    dnames <- dnames[-i]
  }

  attr(ret, "split_vertical") <- split
  attr(ret, "dnames") <- dnames
  
  ## add dimension attributes in ftable-format
  attr(ret, "col.vars") <- dnames[split]
  attr(ret, "row.vars") <- dnames[!split]

  class(ret) <- class(x)
  ret
}

"[<-.structable" <- function(x, ..., value) {
  args <- if (nargs() < 4)
    list(..1)
  else
    list(..1, ..2)
  
  ## handle one-arg cases
  if (nargs() < 4)
    return(if (length(args[[1]]) > 1)
               ## resolve calls like x[c(1,2)]
               Recall(x, args[[1]][1],
                      value = Recall(x[ args[[1]][1] ], args[[1]][-1], value = value))
           else
               ## resolve x[foo] 
               if (attr(x, "split_vertical")[1])
                   Recall(x,,args[[1]], value = value)
               else
                   Recall(x,args[[1]],, value = value)
           )
  
  ## handle calls like x[c(1,2), c(3,4)]
  if (length(args[[1]]) > 1 && length(args[[2]]) > 1)
    return(Recall(x, args[[1]][1], args[[2]][1],
                  value = Recall(x[ args[[1]][1], args[[2]][1] ],
                    args[[1]][-1], args[[2]][-1], value = value)))
  
  ## handle calls like x[c(1,2), 3]
  if (length(args[[1]]) > 1)
    return(Recall(x, args[[1]][1], args[[2]],
                  value = Recall(x[ args[[1]][1], args[[2]] ],
                    args[[1]][-1], ,value = value)))
  
  ## handle calls like x[1, c(1,3)]
  if (length(args[[2]]) > 1)
    return(Recall(x, args[[1]], args[[2]][1],
                  value = Recall(x[ args[[1]], args[[2]][1] ],,
                    args[[2]][-1], value = value)))

  ## final cases like x[1,2] or x[1,] or x[,1]
  dnames <- attr(x, "dnames")
  split <- attr(x, "split_vertical")
  rv <- dnames[!split]
  cv <- dnames[split]

  lsym <- is.symbol(args[[1]])
  rsym <- is.symbol(args[[2]])
  if (!lsym) {
    rstep <- dim(unclass(x))[1] / length(rv[[1]])
    if (is.character(args[[1]]))
      args[[1]] <- match(args[[1]], rv[[1]])
  }
  if (!rsym) {
    cstep <- dim(unclass(x))[2] / length(cv[[1]])
    if (is.character(args[[2]]))
      args[[2]] <- match(args[[2]], cv[[1]])
  }

  lind <- if (!lsym)
    (1 + (args[[1]] - 1) * rstep) : (args[[1]] * rstep)
  else
    1:nrow(unclass(x))
  rind <- if (!rsym)
    (1 + (args[[2]] - 1) * cstep) : (args[[2]] * cstep)
  else
    1:ncol(unclass(x))
  ret <- unclass(x)

  ret[lind, rind] <- value

  class(ret) <- class(x)
  ret
}

as.table.structable <- function(x, ...) {
  ret <- stats:::as.table.ftable(x) 
  structure(aperm(ret, match(names(dimnames(ret)), names(attr(x, "dnames")))),
            class = "table")
}

plot.structable <- function(x, ...)
  mosaic(x, ...)

t.structable <- function(x) {
  ret <- t.default(x)
  attr(ret, "split_vertical") <- !attr(ret, "split_vertical")
  hold <- attr(ret, "row.vars")
  attr(ret, "row.vars") = attr(ret, "col.vars")
  attr(ret, "col.vars") = hold
  ret
}

is.structable <- function(x)
  inherits(x, "structable")

dim.structable <- function(x)
  as.integer(sapply(attr(x, "dnames"), length))

print.structable <- function(x, ...) {
  class(x) <- "ftable"
  NextMethod("print", object = x)
}

dimnames.structable <- function(x) attr(x,"dnames")

as.vector.structable <- function(x, ...)
  as.vector(as.table(x), ...)
  
as.matrix.structable <- function(x) matrix(as.vector(unclass(x)), ncol = attr(x, "dim")[2])
