mosaicpairs <- function(x, ...) UseMethod("mosaicpairs")

"mosaicpairs.formula" <-
function(formula, data = NULL, ...,
         main = deparse(substitute(data)), subset)
{
    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())
    if(inherits(edata, "ftable")
       || inherits(edata, "table")
       || length(dim(edata)) > 2) {
        dat <- as.table(data)
        varnames <- attr(terms(formula), "term.labels")
        if(all(varnames != "."))
            ind <- match(varnames, names(dimnames(dat)))
            if (any(is.na(ind)))
              stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", main))
        mosaicpairs(dat, main = main, ...)
    }
    else {
        if(is.matrix(edata))
            m$data <- as.data.frame(data)
        m$... <- NULL
        m[[1]] <- as.name("model.frame")
        mf <- eval(m, parent.frame())
        mosaicpairs(table(mf), main = main, ...)
    }
}

"mosaicpairs.default" <-
  function(x, main = deparse(substitute(x)), xlab = NULL, ylab = NULL, labels, ...,
           type = c("pairwise", "total", "conditional", "joint"),
           shade = TRUE, oma = NULL, cex.labels = NULL, label.pos = 0.5,
	   font.labels = 1, gap = 1)
  {


    type <- match.arg(type)
    nc<-length(dim(x))
    index<-1:length(dim(x))
    if (nc < 2)
      stop("dimensions less than 2 in the argument to mosaicpairs")

    if (missing(labels)) {
      labels <-  names(dimnames(x))
      if (is.null(labels))
        labels <- paste("var", 1:nc)
    }
    if (is.null(oma)) {
      oma <- c(4, 4, 4, 4)
      if (!is.null(main))
        oma[3] <- 6
    }

    opar <- par(mfrow = c(nc, nc), mar = rep(gap/2, 4),oma=oma)
    on.exit(par(opar))

    for (i in 1:nc)
      for (j in 1:nc) {mfg <- par("mfg")
                       if (i == j) {

                         plot(1,type = "n",axes = FALSE, xlab = "", ylab = "");
                         par(usr = c(0, 1, 0, 1))
                         if (is.null(cex.labels)) {
                           l.wid <- strwidth(labels, "user")
                           cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                         }
                         text(0.5, label.pos, labels[i], cex = cex.labels, font = font.labels);
                         box();
                       }
                       else

                         {
                           switch(type,
                                  pairwise = mosaicplot(margin.table(x, c(j,i)),
                                    main = NULL,
                                    shade=shade, clegend=FALSE, xlab="", ylab="", cex.axis=1 ),
                                  total = mosaicplot(x, shade = shade, clegend = FALSE,
                                    main = NULL,
                                    xlab="", ylab="", cex.axis=1),
                                  conditional = mosaicplot(margin.table(x,
                                    c(j, i, index[!index %in% c(j, i)])), shade=shade,
                                    margin=list(c(j, index[!index %in% c(j, i)]),
                                      c(i, index[!index %in% c(j, i)])),
                                    main = NULL,
                                    clegend=FALSE, xlab="", ylab="", cex.axis=1.2),
                                  joint = mosaicplot(margin.table(x,
                                    c(j, i, index[!index %in% c(j, i)])), shade=shade,
                                    margin=list(c(j, i), c(index[!index %in% c(j, i)])),
                                    main = NULL,
                                    clegend=FALSE, xlab="", ylab="", cex.axis=1.2)
                                  )

                         }
                     }

 if (!is.null(main))
        mtext(main, 3, 3, TRUE, 0.5, cex = 2)
    invisible(NULL)

  }
