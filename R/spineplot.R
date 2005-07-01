spineplot <- function(x, ...) {
  UseMethod("spineplot")
}

spineplot.default <- function(x, split_vertical = TRUE, ...)
{
  mosaic(x, split_vertical = split_vertical, ...)
}

spineplot.formula <- function(formula, data = list(), cut = NULL,  ...)
{
    mf <- model.frame(formula, data = data)
    if(NCOL(mf) != 2) stop("`formula' should specify exactly two variables")
    y <- mf[,1]
    x <- mf[,2]    
    if(!is.factor(y)) stop("dependent variable should be a factor")
        
    if(!is.factor(x)) {
      if(is.null(cut)) cut <- fivenum(x)
      x <- cut(x, cut, include.lowest = TRUE)
    }

    tab <- table(x, y)
    names(dimnames(tab)) <- names(mf)[2:1]
    
    spineplot(tab, ...)
}
