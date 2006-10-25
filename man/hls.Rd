\name{hls}
\alias{hls}
\title{HLS Color Specification}
\description{
  Create a HLS color from specifying hue, luminance and saturation.
}
\usage{
hls(h = 1, l = 0.5, s = 1)
}
\arguments{
  \item{h}{hue value in [0, 1].}
  \item{l}{luminance value in [0, 1].}
  \item{s}{saturation value in [0, 1].}
}
\details{
  HLS colors are a similar specification of colors as HSV
  colors, but using hue/luminance/saturation rather that
  hue/saturation/value.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link{hsv}},
  \code{\link{hcl2hex}},
  \code{\link[colorspace]{polarLUV}}
}
\examples{
## an HLS color wheel
pie(rep(1, 12), col = sapply(1:12/12, function(x) hls(x)))
}
\keyword{hplot}
