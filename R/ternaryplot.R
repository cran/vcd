"ternaryplot" <-
function (x,
          scale = 1,
          dimnames = NULL,
          dimnames.position = c("corner", "edge", "none"),
          dimnames.color = "black",
          id = NULL,
          id.color = "black",
          coordinates = FALSE,
          grid = TRUE,
          grid.color = "gray",
          labels = c("inside", "outside", "none"),
          labels.color = "darkgray",
          border = "black",
          bg = "white",
          pch = 19,
          cex = 1,
          prop.size = FALSE,
          col = "red",
          main = "ternary plot",
          ...)
{
  ## parameter handling
  labels <- match.arg(labels)
  if (grid == TRUE) grid <- "dotted"

  if (coordinates)
    id <- paste("(",round(x[,1] * scale, 1),",",
                    round(x[,2] * scale, 1),",",
                    round(x[,3] * scale, 1),")", sep="")

  dimnames.position <- match.arg(dimnames.position)
  if(is.null(dimnames) && dimnames.position != "none")
    dimnames <- colnames(x)

  if(is.logical(prop.size) && prop.size) prop.size <- 3
  
  ## some error handling
  if(ncol(x) != 3)
    stop("Need a matrix with 3 columns")
  if(any(x) < 0) stop("X must be non-negative")
  s <- rowSums(x)
  if(any(s <= 0)) stop("each row of X must have a positive sum")

  ## rescaling
#  if(max(abs(s - 1)) > 1e-6) {
#    warning("row(s) of X will be rescaled")
    x <- x / s
#  }
  
  ## prepare plot
  top <- sqrt(3) / 2
  par(plt = c(0.06, 0.94, 0.15, 0.87))
  plot.new()
  xlim <- c(-0.03, 1.03)
  ylim <- c(0, top)
  par(usr = c(xlim, ylim),
      oma = c(0, 0, 1, 0)
      )
  plot.window(xlim = xlim, ylim = ylim, asp = 1)
  eps <- 0.01

  ## coordinates of point P(a,b,c): xp = b + c/2, yp = c * sqrt(3)/2

  ## triangle
  polygon(c(0, 0.5, 1), c(0, top, 0), col = bg, xpd = NA, border = border, ...)

  ## title, labeling
  title(main, outer = TRUE, line = -1)
  if (dimnames.position == "corner") {
    axis(1, at = c(-0.03, 1.03), labels = dimnames[1:2], tick = FALSE, font = 2)
    axis(3, at = 0.5, labels = dimnames[3], tick = FALSE, font = 2)
  }
  if (dimnames.position == "edge") {
    shift <- eps * if (labels == "outside") 8 else 0
    text (0.25 - 2 * eps - shift, 0.5 * top + shift, dimnames[2], srt = 60, col = dimnames.color)
    text (0.75 + 3 * eps + shift, 0.5 * top + shift, dimnames[1], srt = -60, col = dimnames.color)
    text (0.5, 0, dimnames[3], pos = 1, offset = 0.5 + 30 * shift, xpd = NA, col = dimnames.color)
  }

  ## grid
  if (is.character(grid))
    for (i in 1:4 * 0.2) {
      ## a - axis
      lines (c(1 - i , (1 - i) / 2), c(0, 1 - i) * top, lty = grid, col = grid.color)
      ## b - axis
      lines (c(1 - i , 1 - i + i / 2), c(0, i) * top, lty = grid, col = grid.color)
      ## c - axis
      lines (c(i/2, 1 - i + i/2), c(i, i) * top, lty = grid, col = grid.color)

      ## grid labels
      if (labels == "inside") {
        text ((1 - i) * 3 / 4 - eps, (1 - i) / 2 * top, i * scale,
              col = labels.color, srt = 120)
        text (1 - i + i / 4 + eps, i / 2 * top - eps, (1 - i) * scale,
              col = labels.color, srt = -120)
        text (0.5, i * top + eps, i * scale, col = labels.color)
      } 
      if (labels == "outside") {
        text ((1 - i) / 2 - 6 * eps, (1 - i) * top, (1 - i) * scale, col = labels.color)
        text (1 - (1 - i) / 2 + 3 * eps, (1 - i) * top + 5 * eps, i * scale,
              srt = -120, col = labels.color)
        text (i + eps, 0, (1 - i) * scale, pos = 1, offset = 1.5,
              srt = 120, xpd = NA, col = labels.color)
      }
    }
  

  ## plot points
  xp <- x[,2] + x[,3] / 2
  yp <- x[,3] * top
  points(xp, yp, pch = pch, col = col,
         cex = if(prop.size) prop.size * (s / max(s)) else cex, ...)

  ## plot 
  if (!is.null(id))
    text (xp, yp, as.character(id), pos = 1, offset = 0.5 * cex, col = id.color)
}


