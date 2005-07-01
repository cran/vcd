\name{rainbow_hcl}
\alias{rainbow_hcl}
\alias{diverge_hcl}
\alias{diverge_hsv}
\title{HCL and HSV Color Palettes}
\description{
  Color palettes based on the HCL and HSV color spaces.
}
\usage{
rainbow_hcl(n, c = 50, l = 70, start = 0, end = 360*(n-1)/n, \dots)

diverge_hcl(n, h = c(260, 0), c = 100, l = c(90, 50), \dots)
diverge_hsv(n, h = c(2/3, 0), s = 1, v = 1, ...)  
}
\arguments{
  \item{n}{the number of colors (\deqn{\ge 1}{>= 1}) to be in the palette.}
  \item{c}{chroma value in the HCL color description.}
  \item{l}{luminance value in the HCL color description.}
  \item{start}{the hue at which the rainbow begins.}
  \item{end}{the hue at which the rainbow ends.}
  \item{h}{hue value in the HCL or HSV color description, has to be
    in [0, 360] for HCL and in [0, 1] for HSV colors.}
  \item{s}{saturation value in the HSV color description.}
  \item{v}{value value in the HSV color description.}
  \item{\dots}{Other arguments passed to \code{\link[colorspace]{hcl}}
    or \code{\link{hsv}}, respectively.}
}
\details{
  All three functions compute palettes based on either the HCL or
  the HSV color space. \code{rainbow_hcl} computes a rainbow of colors
  defined by different hues given a single value of each chroma and luminance.
  It corresponds to \code{\link{rainbow}} which computes a rainbow in
  HSV space.
  
  Both functions, \code{diverge_hcl} and \code{diverge_hsv}, compute 
  a set of colors diverging from a neutral center (grey or white, without
  color) to two different extreme colors (blue and red by default). This
  is similar to \code{\link{cm.colors}}. For the diverging HSV colors,
  two hues \code{h} are needed, a maximal saturation \code{s} and a
  fixed value \code{v}. The saturation is then varied to obtain the
  diverging colors. For the diverging HCL colors, again two hues \code{h}
  are needed, a maximal chroma \code{chroma} and two luminances \code{l}.
  The colors are then created by linear interpolation between the full
  chroma and minimal luminance and a chroma of 0 and maximal luminance.

  The diverging palettes are used for choosing the colors for mosaic 
  displays in \code{\link{mosaic}}.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link[colorspace]{hcl}},
  \code{\link{hsv}},
  \code{\link{shadings}}}
}
\examples{
pie(rep(1, 12), col = rainbow_hcl(12))
pie(rep(1, 5), col = diverge_hcl(5))
pie(rep(1, 5), col = diverge_hsv(5))
}
\keyword{hplot}
