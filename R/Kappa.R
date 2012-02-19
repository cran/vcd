Kappa <- function (x, weights = c("Equal-Spacing", "Fleiss-Cohen"))
{
  if (is.character(weights))
      weights <- match.arg(weights)

  d  <- diag(x)
  n  <- sum(x)
  nc <- ncol(x)
  colFreqs <- colSums(x)/n
  rowFreqs <- rowSums(x)/n
  
  ## Kappa
  kappa <- function (po, pc)
    (po - pc) / (1 - pc)
  std  <- function (po, pc, W = 1)
    sqrt(sum(W * W * po * (1 - po)) / crossprod(1 - pc) / n)
    
  ## unweighted
  po <- sum(d) / n
  pc <- crossprod(colFreqs, rowFreqs)
  k <- kappa(po, pc)
  s <- std(po, pc)
  
  ## weighted 
  W <- if (is.matrix(weights))
    weights
  else if (weights == "Equal-Spacing")
    1 - abs(outer(1:nc, 1:nc, "-")) / (nc - 1)
  else
    1 - (abs(outer(1:nc, 1:nc, "-")) / (nc - 1))^2
  pow <- sum(W * x) / n
  pcw <- sum(W * colFreqs %o% rowFreqs)
  kw <- kappa(pow, pcw)
  sw <- std(x / n, 1 - pcw, W)

  structure(
            list(Unweighted = c(
                   value = k,
                   ASE   = s
                   ),
                 Weighted = c(
                   value = kw,
                   ASE   = sw
                   ),
                 Weights = W
                 ),
            class = "Kappa"
       )
}

print.Kappa <- function (x, ...) {
  tab <- rbind(x$Unweighted, x$Weighted)
  rownames(tab) <- names(x)[1:2]
  print(tab, ...)
  invisible(x)
}

summary.Kappa <- function (object, ...)
  structure(object, class = "summary.Kappa")

print.summary.Kappa <- function (x, ...) {
  print.Kappa(x, ...)
  cat("\nWeights:\n")
  print(x$Weights, ...)
  invisible(x)
}

confint.Kappa <- function(object, parm, level = 0.95, ...) {
  q <- qnorm((1 + level) / 2)
  matrix(c(object[[1]][1] - object[[1]][2] * q,
           object[[1]][1] + object[[1]][2] * q,
           object[[2]][1] - object[[2]][2] * q,
           object[[2]][1] + object[[2]][2] * q),
         ncol = 2, byrow = TRUE, 
         dimnames = list(Kappa = c("Unweighted","Weighted"), c("lwr","upr"))
         )
}

