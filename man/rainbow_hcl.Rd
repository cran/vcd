\name{rainbow_hcl}
\alias{rainbow_hcl}
\alias{diverge_hcl}
\alias{diverge_hsv}
\alias{heat_hcl}
\alias{sequential_hcl}
\alias{terrain_hcl}
\alias{hcl2hex}
\title{HCL and HSV Color Palettes}
\encoding{latin1}
\description{
  Color palettes based on the HCL and HSV color spaces.
}
\usage{
rainbow_hcl(n, c = 50, l = 70, start = 0, end = 360*(n-1)/n, \dots)

sequential_hcl(n, h = 260, c. = c(80, 0), l = c(30, 90), power = 1.5, \dots)
heat_hcl(n, h = c(0, 90), c. = c(100, 30), l = c(50, 90), power = c(1/5, 1), \dots)
terrain_hcl(n, h = c(130, 0), c. = c(80, 0), l = c(60, 95), power = c(1/10, 1), \dots)

diverge_hcl(n, h = c(260, 0), c = 80, l = c(30, 90), power = 1.5, \dots)
diverge_hsv(n, h = c(2/3, 0), s = 1, v = 1, power = 1, \dots)

hcl2hex(h = 0, c = 35, l = 85, gamma = 2.2, fixup = TRUE)
}
\arguments{
  \item{n}{the number of colors (\eqn{\ge 1}{>= 1}) to be in the palette.}
  \item{c, c.}{chroma value in the HCL color description.}
  \item{l}{luminance value in the HCL color description.}
  \item{start}{the hue at which the rainbow begins.}
  \item{end}{the hue at which the rainbow ends.}
  \item{h}{hue value in the HCL or HSV color description, has to be
    in [0, 360] for HCL and in [0, 1] for HSV colors.}
  \item{s}{saturation value in the HSV color description.}
  \item{v}{value value in the HSV color description.}
  \item{power}{control parameter determining how chroma
    and luminance should be increased (1 = linear, 2 = quadratic, etc.).}
  \item{gamma}{gamma value of the display.}
  \item{fixup}{logical. Should the color be corrected to a valid RGB
    value before correction?}
  \item{\dots}{Other arguments passed to \code{hcl2hex} or \code{\link{hsv}}, respectively.}
}
\details{
  All functions compute palettes based on either the HCL or
  the HSV color space. For the latter, the implementation \code{\link{hsv}}
  from base R is used. The former is implemented in the \pkg{colorspace}
  package: HCL color space is obtained by taking polar coordinates
  in CIELUV space which is available via \code{\link[colorspace]{polarLUV}}
  from the \pkg{colorspace} package. To transform the HCL coordinates to
  a hexadecimal color string (as returned by \code{hsv}), the function
  \code{\link[colorspace]{hex}} is employed. A convenience wrapper \code{hcl2hex}
  is provided.
  
  \code{rainbow_hcl} computes a rainbow of colors
  defined by different hues given a single value of each chroma and luminance.
  It corresponds to \code{\link{rainbow}} which computes a rainbow in
  HSV space.

  \code{sequential_hcl} gives a sequential palette starting at the full
  color HCL(\code{h}, \code{c[1]}, \code{l[1]}) through to a light color
  HCL(\code{h}, \code{c[2]}, \code{l[2]}) by interpolation.
  
  \code{diverge_hcl} and \code{diverge_hsv}, compute 
  a set of colors diverging from a neutral center (grey or white, without
  color) to two different extreme colors (blue and red by default). This
  is similar to \code{\link{cm.colors}}. For the diverging HSV colors,
  two hues \code{h} are needed, a maximal saturation \code{s} and a
  fixed value \code{v}. The saturation is then varied to obtain the
  diverging colors. For the diverging HCL colors, again two hues \code{h}
  are needed, a maximal chroma \code{c} and two luminances \code{l}.
  The colors are then created by an interpolation between the 
  full color HCL(\code{h[1]}, \code{c}, \code{l[1]}), a neutral color HCL(\code{h}, 0, \code{l[2]})
  and the other full color HCL(\code{h[2]}, \code{c}, \code{l[1]}).

  The palette \code{heat_hcl} gives an implementation of \code{\link{heat.colors}}
  in HCL space. By default, it goes from a red to a yellow hue, while
  simultaneously going to lighter colors (i.e., increasing luminance)
  and reducing the amount of color (i.e., decreasing chroma).
  The \code{terrain_hcl} palette simply calls \code{heat_hcl} with different
  parameters, providing colors similar in spirit to \code{terrain.colors}.
  The lighter colors are not strictly HCL colors, though.

  The diverging palettes are used for choosing the colors for mosaic 
  displays in \code{\link{mosaic}}.
}

\references{
  Zeileis A. and Hornik K. (2006), Choosing Color Palettes for Statistical Graphics.
  Report 41, Department of Statistics and Mathematics, Wirtschaftsuniversität Wien,
  Research Report Series, \url{http://epub.wu-wien.ac.at/}.
}

\value{
  A character vector with RGB codings of the colors in the palette.
}

\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link[colorspace]{polarLUV}},
  \code{\link[colorspace]{hex}},
  \code{\link{hsv}},
  \code{\link{shadings}}}
}
\examples{
## convenience demo functions
wheel <- function(col, radius = 1, ...)
  pie(rep(1, length(col)), col = col, radius = radius, ...) 

pal <- function(col, border = "light gray")
{
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab="")
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

## qualitative palette
wheel(rainbow_hcl(12))

## a few useful diverging HCL palettes
par(mar = rep(0, 4), mfrow = c(4, 1))
pal(diverge_hcl(7))
pal(diverge_hcl(7, h = c(246, 40), c = 96, l = c(65, 90)))
pal(diverge_hcl(7, h = c(130, 43), c = 100, l = c(70, 90)))
pal(diverge_hcl(7, h = c(180, 70), c = 70, l = c(90, 95)))
pal(diverge_hcl(7, h = c(180, 330), c = 59, l = c(75, 95)))
pal(diverge_hcl(7, h = c(128, 330), c = 98, l = c(65, 90)))
pal(diverge_hcl(7, h = c(255, 330), l = c(40, 90)))
pal(diverge_hcl(7, c = 100, l = c(50, 90), power = 1))

## sequential palettes
pal(sequential_hcl(12))
pal(heat_hcl(12, h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1))
pal(terrain_hcl(12, c = c(65, 0), l = c(45, 95), power = c(1/3, 1.5)))
pal(heat_hcl(12, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))

## compare base and vcd palettes
par(mfrow = c(2, 1))
wheel(rainbow(12));wheel(rainbow_hcl(12))
pal(diverge_hcl(7, c = 100, l = c(50, 90)));pal(diverge_hsv(7))
pal(diverge_hcl(7, h = c(180, 330), c = 59, l = c(75, 95)));pal(cm.colors(7))
pal(heat_hcl(12));pal(heat.colors(12))
pal(terrain_hcl(12));pal(terrain.colors(12))
}
\keyword{hplot}
