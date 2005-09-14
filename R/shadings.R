## shading-generating functions should take at least the arguments
##   observed, residuals, expected, df
## and return a function which takes a single argument (interpreted
## to be a vector of residuals).

shading_hsv <- function(observed, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), s = c(1, 0), v = c(1, 0.5),
  interpolate = c(2, 4), lty = 1, eps = NULL, line_col = "black",
  p.value = NULL, level = 0.95, ...)
{
  ## get h/s/v and lty
  my.h <- rep(h, length.out = 2)  ## positive and negative hue
  my.s <- rep(s, length.out = 2)  ## maximum and minimum saturation
  my.v <- rep(v, length.out = 2)  ## significant and non-significant value
  lty <- rep(lty, length.out = 2) ## positive and negative lty

  ## model fitting (if necessary)
  if(is.null(expected) && !is.null(residuals)) stop("residuals without expected values specified")
  if(!is.null(expected) && is.null(df) && is.null(p.value)) {
    warning("no default inference available without degrees of freedom")
    p.value <- NA
  }
  if(is.null(expected) && !is.null(observed)) {
    expected <- loglin(observed, 1:length(dim(observed)), fit = TRUE, print = FALSE)
    df <- expected$df
    expected <- expected$fit
  }
  if(is.null(residuals) && !is.null(observed)) residuals <- (observed - expected)/sqrt(expected)
    
  ## conduct significance test (if specified)
  if(is.null(p.value)) p.value <- function(observed, residuals, expected, df)
    pchisq(sum(as.vector(residuals)^2), df, lower.tail = FALSE)
  if(!is.function(p.value) && is.na(p.value)) {
    v <- my.v[1]
    p.value <- NULL
  } else {
    if(is.function(p.value)) p.value <- p.value(observed, residuals, expected, df)
    v <- if(p.value < (1-level)) my.v[1] else my.v[2]
  }

  ## set up function for interpolation of saturation
  if(!is.function(interpolate)) {
    col.bins <- sort(interpolate)
    interpolate <- stepfun(col.bins,  seq(my.s[2], my.s[1], length = length(col.bins) + 1))
    col.bins <- sort(unique(c(col.bins, 0, -col.bins)))
  } else {
    col.bins <- NULL
  }

  ## store color and lty information for legend
  legend <- NULL
  if(!is.null(col.bins)) {
    res2 <- col.bins
    res2 <- c(head(res2, 1) - 1, res2[-1] - diff(res2)/2, tail(res2, 1) + 1)
    legend.col <- hsv(ifelse(res2 > 0, my.h[1], my.h[2]),
                      pmax(pmin(interpolate(abs(res2)), 1), 0),
		      v, ...)
    lty.bins <- 0
    legend.lty <- lty[2:1]
    legend <- list(col = legend.col, col.bins = col.bins,
                   lty = legend.lty, lty.bins = lty.bins)
  }

  ## set up function that computes color/lty from residuals
  rval <- function(x) {
    res <- as.vector(x)

    fill <- hsv(ifelse(res > 0, my.h[1], my.h[2]),
                pmax(pmin(interpolate(abs(res)), 1), 0),
	        v, ...)
    dim(fill) <- dim(x)
    
    col <- rep(line_col, length.out = length(res))
    if(!is.null(eps)) {
      eps <- abs(eps)
      col[res > eps] <- hsv(my.h[1], 1, v, ...)
      col[res < -eps] <- hsv(my.h[2], 1, v, ...)      
    }
    dim(col) <- dim(x)
    
    lty <- ifelse(x > 0, lty[1], lty[2])    
    dim(lty) <- dim(x)

    return(structure(list(col = col, fill = fill, lty = lty), class = "gpar"))
  }
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- p.value
  return(rval)
}
class(shading_hsv) <- "grapcon_generator"


shading_hcl <- function(observed, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL,
  interpolate = c(2, 4), lty = 1, eps = NULL, line_col = "black",
  p.value = NULL, level = 0.95, ...)
{
  ## set defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- c(100, 20)
  if(is.null(l)) l <- c(90, 50)

  ## get h/c/l and lty
  my.h <- rep(h, length.out = 2)  ## positive and negative hue
  my.c <- rep(c, length.out = 2)  ## significant and non-significant maximum chroma
  my.l <- rep(l, length.out = 2)  ## maximum and minimum luminance
  lty <- rep(lty, length.out = 2) ## positive and negative lty

  ## model fitting (if necessary)
  if(is.null(expected) && !is.null(residuals)) stop("residuals without expected values specified")
  if(!is.null(expected) && is.null(df) && is.null(p.value)) {
    warning("no default inference available without degrees of freedom")
    p.value <- NA
  }
  if(is.null(expected) && !is.null(observed)) {
    expected <- loglin(observed, 1:length(dim(observed)), fit = TRUE, print = FALSE)
    df <- expected$df
    expected <- expected$fit
  }
  if(is.null(residuals) && !is.null(observed)) residuals <- (observed - expected)/sqrt(expected)
    
  ## conduct significance test (if specified)
  if(is.null(p.value)) p.value <- function(observed, residuals, expected, df)
    pchisq(sum(as.vector(residuals)^2), df, lower.tail = FALSE)
  if(!is.function(p.value) && is.na(p.value)) {
    max.c <- my.c[1]
    p.value <- NULL
  } else {
    if(is.function(p.value)) p.value <- p.value(observed, residuals, expected, df)
    max.c <- ifelse(p.value < (1-level), my.c[1], my.c[2])
  }

  ## set up function for interpolation of saturation
  if(!is.function(interpolate)) {
    col.bins <- sort(interpolate)
    interpolate <- stepfun(col.bins,  seq(0, 1, length = length(col.bins) + 1))
    col.bins <- sort(unique(c(col.bins, 0, -col.bins)))
  } else {
    col.bins <- NULL
  }

  ## store color and lty information for legend
  legend <- NULL
  if(!is.null(col.bins)) {
    res2 <- col.bins
    res2 <- c(head(res2, 1) - 1, res2[-1] - diff(res2)/2, tail(res2, 1) + 1)
    legend.col <- hcl(ifelse(res2 > 0, my.h[1], my.h[2]),
                      max.c * pmax(pmin(interpolate(abs(res2)), 1), 0),
	              my.l[1] + diff(my.l) * pmax(pmin(interpolate(abs(res2)), 1), 0),
		      ...)
    lty.bins <- 0
    legend.lty <- lty[2:1]
    legend <- list(col = legend.col, col.bins = col.bins,
                   lty = legend.lty, lty.bins = lty.bins)
  }

  ## set up function that computes color/lty from residuals
  rval <- function(x) {
    res <- as.vector(x)

    fill <- hcl(ifelse(res > 0, my.h[1], my.h[2]),
                max.c * pmax(pmin(interpolate(abs(res)), 1), 0),
	        my.l[1] + diff(my.l) * pmax(pmin(interpolate(abs(res)), 1), 0),
	        ...)
    dim(fill) <- dim(x)
    
    col <- rep(line_col, length.out = length(res))
    if(!is.null(eps)) {
      eps <- abs(eps)
      col[res > eps] <- hcl(my.h[1], max.c, my.l[2], ...)
      col[res < -eps] <- hcl(my.h[2], max.c, my.l[2], ...)
    }
    dim(col) <- dim(x)
    
    lty <- ifelse(x > 0, lty[1], lty[2])    
    dim(lty) <- dim(x)

    return(structure(list(col = col, fill = fill, lty = lty), class = "gpar"))
  }
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- p.value
  return(rval)
}
class(shading_hcl) <- "grapcon_generator"

shading_Friendly <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), lty = 1:2, interpolate = c(2, 4), eps = 0.01, line_col = "black", ...)
{
  shading_hsv(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
              h = h, v = 1, lty = lty, interpolate = interpolate,
	      eps = eps, line_col = line_col, p.value = NA, ...)
}
class(shading_Friendly) <- "grapcon_generator"

shading_max <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL, lty = 1, eps = NULL, line_col = "black", level = c(0.9, 0.99), n = 1000, ...)
{
  stopifnot(length(dim(observed)) == 2)
  
  ## set defaults
  if(is.null(h)) h <- c(260, 0)
  if(is.null(c)) c <- c(100, 20)
  if(is.null(l)) l <- c(90, 50)  
  
  obs.test <- coindep_test(observed, n = n)
  col.bins <- obs.test$qdist(sort(level))
  rval <- shading_hcl(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
                        h = h, c = c, l = l, interpolate = col.bins, lty = lty,
			eps = eps, line_col = line_col, p.value = obs.test$p.value, ...)
  return(rval)
}
class(shading_max) <- "grapcon_generator"

shading_binary <- function(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  col = hcl(c(260, 0), 50, 70))
{
  ## check col argument
  if (length(col) != 2) stop("Need exactly two colors!")
  
  ## store color information for legend
  legend <- list(col = col[2:1], col.bins = 0, lty = NULL, lty.bins = NULL)

  ## set up function that computes color/lty from residuals
  rval <- function(x)
    gpar(fill = ifelse(x > 0, col[1], col[2]))

  ## add meta information for legend
  attr(rval, "legend") <- legend
  attr(rval, "p.value") <- NULL
  
  rval
}
class(shading_binary) <- "grapcon_generator"



## color palettes

rainbow_hcl <- function(n, c = 50, l = 70, start = 0, end = 360*(n-1)/n, ...)
{
  if(n > 0) hcl(seq(start, end, length = n), c = c, l = l, ...)
    else character(0)
}

diverge_hcl <- function(n, h = c(260, 0), c = 100, l = c(90, 50), ...)
{
  if(n < 1) return(character(0))
  h <- rep(h, length.out = 2)
  c <- c[1]
  l <- rep(l, length.out = 2)
  rval <- seq(1, -1, length = n)
  rval <- hcl(h = ifelse(rval > 0, h[1], h[2]),
              c = c * abs(rval),
              l = l[1] + diff(l) * abs(rval),
              ...)
  return(rval)
}

diverge_hsv <- function(n, h = c(2/3, 0), s = 1, v = 1, ...)
{
  if(n < 1) return(character(0))
  h <- rep(h, length.out = 2)
  s <- s[1]
  v <- v[1]
  rval <- seq(-s, s, length = n)
  rval <- hsv(h = ifelse(rval > 0, h[2], h[1]), s = abs(rval), v = v, ...)
  return(rval)
}
