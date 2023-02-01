odds <- function(x, log = FALSE, ...)
	lodds(x, log = log, ...)


lodds <- function(x, ...)
    UseMethod("lodds")

lodds.formula <-
    function(formula, data = NULL, ..., subset = NULL, na.action = NULL)
{
    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
    vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
    varnames <- vars[[1]]

    condnames <- if (length(vars) > 1) vars[[2]] else NULL

    dep <- gsub(" ", "", fstr[[1]][1])
    if (!dep %in% c("","Freq")) {
        if (all(varnames == ".")) {
            varnames <- if (is.data.frame(data))
                            colnames(data)
                        else
                            names(dimnames(as.table(data)))
            varnames <- varnames[-which(varnames %in% dep)]
        }

        varnames <- c(dep, varnames)
    }


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
                ind <- c(ind, condind)
            }
            dat <- margin.table(dat, ind)
        }
        lodds.default(dat, strata = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
    } else {
        m <- m[c(1, match(c("formula", "data", "subset", "na.action"), names(m), 0))]
        m[[1]] <- as.name("xtabs")
        m$formula <-
            formula(paste(if("Freq" %in% colnames(data)) "Freq",
                          "~",
                          paste(c(varnames, condnames), collapse = "+")))
        tab <- eval(m, parent.frame())
        lodds.default(tab, ...)
    }
}

lodds.default <- function(x, response = NULL, strata = NULL, log = TRUE,
                                  ref = NULL, correct = any(x == 0), ...)
{
    ## check dimensions
    L <- length(d <- dim(x))
    if(any(d < 2L)) stop("All table dimensions must be 2 or greater")

		## assign and check response and stata; convert variable names to indices
    if (is.null(response)) {
    	if (is.null(strata)) {
    		response <- 1
    		strata <- setdiff(1:L, response)
    	}
    	else {   # only strata was specified
    		if(L - length(strata) != 1L) stop("All but 1 dimension must be specified as strata.")
    		if(is.character(strata)) strata <- which(names(dimnames(x)) == strata)
    		response <- setdiff(1:L, strata)
    	}
    }
    else {  # response was specified; take strata as the complement
    	if(length(response) > 1) stop("Only 1 dimension can be specified as a response")
    	if(is.character(response)) response <- which(names(dimnames(x)) == response)
    	if (!is.null(strata))
	    	warning(paste("strata =", paste(strata, collapse=","), "ignored when response has been specified"))
    	strata <- setdiff(1:L, response)
	  }

    ## dimensions of primary R x 1 table  ### [Or should this just be a vector???]
    dp <- if (length(strata)) d[response] else d
    dn <- if (length(strata)) dimnames(x)[response] else dimnames(x)
    R <- dp[1]
    C <- 1
                                        # shadow matrix with proper dimnames
    X <- matrix(0, R, C, dimnames=dn)

    ## process reference category
    if(!is.null(ref)) {
        if(is.character(ref)) {
            ref <- match(ref, colnames(x))
        } else if(is.numeric(ref)) {
            ref <- as.integer(ref)
        } else {
            stop("Wrong 'ref=' argument!")
        }
    }

    ## compute corresponding indices
    compute_index <- function(n, ref) {
        if(is.null(ref)) return(cbind(1:(n-1), 2:n))
        rval <- cbind(ref, 1:n)
        d <- rval[,2L] - rval[,1L]
        rval <- rbind(
            rval[d > 0, 1:2],
            rval[d < 0, 2:1]
        )
        return(rval[order(rval[,1L]),,drop = FALSE])
    }
    Rix <- compute_index(R, ref[[1L]])
    contr <- matrix(0L, nrow = (R-1), ncol = R)
    colnames(contr) <- rownames(X)
    rownames(contr) <- rep("", R-1)
    for(i in 1:(R-1)) {
        rix <-  i
        cix <- Rix[i,]
        contr[rix, cix] <- c(1L, -1L)
        rownames(contr)[rix] <-  paste(rownames(X)[Rix[i,]], collapse = ":")
    }

    ## handle strata
    if (!is.null(strata)) {
  	if (length(strata)==1) {
            sn <- dimnames(x)[[strata]]
        }
        else {
            sn <- apply(expand.grid(dimnames(x)[strata]), 1, paste, collapse = ":")
        }
        rn <- as.vector(outer( dimnames(contr)[[1]], sn, paste, sep='|'))
        cn <- as.vector(outer( dimnames(contr)[[2]], sn, paste, sep='|'))
        contr <- kronecker(diag(prod(dim(x)[strata])), contr)
        rownames(contr) <- rn
        colnames(contr) <- cn
    }

    ## dimnames for array version
    dn <- list(rep("", R-1))
    for(i in 1:(R-1)) dn[[1]][i] <- paste(rownames(X)[Rix[i,]], collapse = ":")
    if (!is.null(strata)) dn <- c(dn, dimnames(x)[strata])
    ndn <- names(dimnames(x))
    if (!is.null(names(dimnames(x)))) names(dn) <- c(ndn[response], ndn[strata])

    ## point estimates

    if (is.logical(correct)) {
        add <- if(correct) 0.5 else 0
    }
    else if(is.numeric(correct)) {
        add <- as.vector(correct)
        if (length(add) != length(x))
            stop("array size of 'correct' does not conform to the data")
    }
    else stop("correct is not valid")

    ## reorder columns of contrast matrix to match original data
    contr <- contr[, order(as.vector(aperm(array(seq.int(prod(d)), d), c(response, strata))))]

    ##coef <- drop(contr %*% log(as.vector(x) + add))
    ##FIXME: 0 cells mess up the matrix product, try workaround:
    mat <- log(as.vector(x) + add) * t(contr)
    nas <- apply(contr != 0 & is.na(t(mat)), 1, any)
    coef <- apply(mat, 2, sum, na.rm = TRUE)
    coef[nas] <- NA
    ## covariances
    ##vcov <- crossprod(diag(sqrt(1/(as.vector(x) + add))) %*% t(contr))
    tmp <- sqrt(1/(as.vector(x) + add)) * t(contr)
    tmp[is.na(tmp)] <- 0
    vcov <- crossprod(tmp)
    vcov[nas,] <- NA
    vcov[,nas] <- NA

    rval <- structure(list(
        response = response,
        strata = strata,
        coefficients = coef,
        dimnames = dn,
        dim = as.integer(sapply(dn, length)),
        vcov = vcov,
        contrasts = contr,
        log = log
    ), class = "lodds")
    rval
}

## ----------------  Methods -------------------

summary.lodds <- function(object, ...)
    lmtest::coeftest(object, ...)

## dim methods
dimnames.lodds <- function(x, ...) x$dimnames
dim.lodds <- function(x, ...) x$dim

## t/aperm-methods
t.lodds <- function(x)
    aperm(x)

### FIXME
aperm.lodds <- function(a, perm = NULL, ...)
{
    d <- length(a$dim)
    if(is.null(perm)) {
        perm <- if (d < 3) 2L : 1L else c(2L : 1L, d : 3L)
    } else {
        if (any(perm[1:2] > 2L) || (d > 2L) && any(perm[-c(1:2)] < 2L))
            stop("Mixing of strata and non-strata variables not allowed!")
    }
    nams <- names(a$coefficients)
    a$coefficients <- as.vector(aperm(array(a$coef, dim = a$dim),
                                      perm, ...))
    nams <- as.vector(aperm(array(nams, dim = a$dim), perm, ...))
    names(a$coefficients) <- nams
    a$dimnames <- a$dimnames[perm]
    a$dim <- a$dim[perm]
    a$vcov <- a$vcov[nams, nams]
    a$contrasts <- a$contrasts[nams,]
    a
}


## straightforward methods
coef.lodds <- function(object, log = object$log, ...)
    if(log) object$coefficients else exp(object$coefficients)

vcov.lodds <- function(object, log = object$log, ...)
    if(log) object$vcov else `diag<-`(object$vcov, diag(object$vcov) * exp(object$coefficients)^2)

confint.lodds <-
    function(object, parm, level = 0.95, log = object$log, ...) {
        if (log) confint.default(object, parm = parm, level = level, ... )
        else {
            object$log = TRUE
            exp(confint.default(object, parm = parm, level = level, ... ))
        }
    }


### DONE:
##  The header should say:
#    (log) odds for vn[response] by ... all the rest (strata)
#   Fixed:  clash with make_header in loddsratio
make_header_odds <- function(x)
{
    vn <- names(dimnames(x))
    resp <- vn[x$response]
    strat <- paste(vn[x$strata], collapse=", ")
    header <- c(if(x$log) "log" else "",
                "odds for", resp, "by", strat,
#                if (length(vn)>2) c("by", paste(vn[-(1:2)], collapse=', ')),
                "\n\n")
    paste(header, sep = " ")
}

## print method
print.lodds <- function(x, log = x$log, ...) {
    cat(make_header_odds(x))
    print(drop(array(coef(x, log = log), dim = dim(x), dimnames = dimnames(x)), ...))
    invisible(x)
}

## as.data.frame
#as.data.frame.lodds <-
#    function(x, ...)
#        as.data.frame.table(vcd:::as.array.loddsratio(x), ...)

## Q:  I don't understand the purpose of the row.names and optional arguments
## DM: The generic has them, so each method must have them, too
as.data.frame.lodds <- function(x, row.names = NULL, optional = FALSE, log=x$log, ...) {
    df <-data.frame(expand.grid(dimnames(x)),
                    logodds = coef(x, log = log),
                    ASE = sqrt(diag(vcov(x, log = log))), row.names = row.names, ...
                    )
    if (!log) colnames(df)[ncol(df) - 1] <- "odds"
    df
}

## FIXME
## reshape coef() methods
as.matrix.lodds <- function (x, log=x$log, ...) {
    ## Coef <- coef(x, log = log)
    ## if (length(dim(x))==2) matrix(Coef, ncol = dim(x)[2], dimnames=dimnames(x))
    ## else {  # drop leading dimensions with length 1, then reshape
    ##     ddim <- which(dim(x)[1:2]==1)
    ##     dim(Coef) <- dim(x)[-ddim]
    ##     dimnames(Coef) <- dimnames(x)[-ddim]
    ##     if (length(dim(Coef))==1) Coef
    ##     else
    ##         matrix(Coef, ncol = prod(dim(Coef)[-1]),
    ##                dimnames=list(dimnames(Coef)[[1]], apply(expand.grid(dimnames(Coef)[[-1]]), 1, paste, collapse = ":")))
    ## }
    as.array(x, log = log, ...)
}


as.array.lodds <- function (x, log=x$log, ...) {
    res <- array(coef(x, log = log), dim = dim(x), dimnames=dimnames(x))
    drop(res)
}
