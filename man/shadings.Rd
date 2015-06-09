\name{shadings}
\alias{shadings}
\alias{shading_hsv}
\alias{shading_hcl}
\alias{shading_max}
\alias{shading_Friendly}
\alias{shading_Friendly2}
\alias{shading_Marimekko}
\alias{shading_diagonal}
\alias{shading_sieve}
\alias{shading_binary}
\alias{hcl2hex}
\encoding{UTF-8}
\title{Shading-generating Functions for Residual-based Shadings}
\description{
  Shading-generating functions for computing residual-based shadings 
  for mosaic and association plots.
}
\usage{
shading_hcl(observed, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL, interpolate = c(2, 4), lty = 1,
  eps = NULL, line_col = "black", p.value = NULL, level = 0.95, \dots)

shading_hsv(observed, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), s = c(1, 0), v = c(1, 0.5),
  interpolate = c(2, 4), lty = 1, eps = NULL, line_col = "black",
  p.value = NULL, level = 0.95, \dots)

shading_max(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = NULL, c = NULL, l = NULL, lty = 1, eps = NULL, line_col = "black",
  level = c(0.9, 0.99), n = 1000, \dots)

shading_Friendly(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = c(2/3, 0), lty = 1:2, interpolate = c(2, 4),
  eps = 0.01, line_col = "black", \dots)

shading_Friendly2(observed = NULL, residuals = NULL,
  expected = NULL, df = NULL, lty = 1:2, interpolate = c(2, 4), eps =
  0.01, line_col = "black", \dots)

shading_sieve(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  h = c(260, 0), lty = 1:2, interpolate = c(2, 4),
  eps = 0.01, line_col = "black", \dots)

shading_binary(observed = NULL, residuals = NULL, expected = NULL, df = NULL,
  col = NULL)

shading_Marimekko(x, fill = NULL, byrow = FALSE)

shading_diagonal(x, fill = NULL)

hcl2hex(h = 0, c = 35, l = 85, fixup = TRUE)
}
\arguments{
  \item{observed}{contingency table of observed values}
  \item{residuals}{contingency table of residuals}
  \item{expected}{contingency table of expected values}
  \item{df}{degrees of freedom of the associated independence model.}  
  \item{h}{hue value in the HCL or HSV color description, has to be
    in [0, 360] for HCL and in [0, 1] for HSV colors. The default is
    to use blue and red for positive and negative residuals respectively.
    In the HCL specification it is \code{c(260, 0)} by default and for HSV
    \code{c(2/3, 0)}.}
  \item{c}{chroma value in the HCL color description. This controls the maximum
    chroma for significant and non-significant results respectively and defaults 
    to \code{c(100, 20)}.}
  \item{l}{luminance value in the HCL color description. Defaults to \code{c(90, 50)}
    for small and large residuals respectively.}
  \item{s}{saturation value in the HSV color description. Defaults to \code{c(1, 0)}
    for large and small residuals respectively.}
  \item{v}{saturation value in the HSV color description. Defaults to \code{c(1, 0.5)}
    for significant and non-significant results respectively.}
  \item{interpolate}{a specification for mapping the absolute size of the residuals to
    a value in [0, 1]. This can be either a function or a numeric vector. In the latter
    case, a step function with steps of equal size going from 0 to 1 is used.}
  \item{lty}{a vector of two line types for positive and negative residuals respectively.
    Recycled if necessary.}
  \item{eps}{numeric tolerance value below which absolute residuals are considered to be
    zero, which is used for coding the border color and line type. If set to \code{NULL} (default),
    all borders have the default color specified by \code{line\_col}.
    If set to a numeric value, all border colors corresponding
    to residuals with a larger absolute value are set to the full positive or negative
    color, respectively;  borders corresponding to smaller residuals are are drawn with
    \code{line\_col} and \code{lty[1]}}.  This is used principally in \code{shading\_Friendly}.
  \item{line_col}{default border color (for \code{shading_sieve}: default sieve color).}
  \item{p.value}{the \eqn{p} value associated with the independence model. By default,
    this is computed from a Chi-squared distribution with \code{df} degrees of freedom.
    \code{p.value} can be either a scalar or a \code{function(observed, residuals, expected, df)}
    that computes the \eqn{p} value from the data. If set to \code{NA} no inference
    is performed.}
  \item{level}{confidence level of the test used. If \code{p.value} is smaller than
    \code{1 - level}, bright colors are used, otherwise dark colors are employed. For
    \code{shading_max} a vector of levels can be supplied. The corresponding critical
    values are then used as \code{interpolate} cut-offs.}
  \item{n}{number of permutations used in the call to \code{coindep_test}.}
  \item{col}{a vector of two colors for positive and negative residuals respectively.}
  \item{fixup}{logical. Should the color be corrected to a valid RGB
    value before correction?}
  \item{x}{object of class \code{table} used to determine the
    dimension.}
  \item{fill}{Either a character vector of color codes, or a palette
    function that generates such a vector. Defaults to
    \code{\link[colorspace]{rainbow_hcl}}}
  \item{byrow}{logical; shall tiles be filled by row or by column?}
  \item{\dots}{Other arguments passed to \code{\link{hcl2hex}}
    or \code{\link{hsv}}, respectively.}
}
\details{
  These shading-generating functions can be passed to \code{strucplot} to generate
  residual-based shadings for contingency tables. \code{strucplot} calls these
  functions with the arguments \code{observed}, \code{residuals}, \code{expected},
  \code{df} which give the observed values, residuals, expected values and associated
  degrees of freedom for a particular contingency table and associated independence 
  model.
  
  The shadings \code{shading_hcl} and \code{shading_hsv} do the same thing conceptually,
  but use HCL or HSV colors respectively. The former is usually preferred because they
  are perceptually based. Both shadings visualize the \emph{sign} of the residuals of
  an independence model using two hues (by default: blue and red). The \emph{absolute size} of 
  the residuals is visualized by the colorfulness and the amount of grey, by default in three categories:
  very colorful for large residuals (> 4), less colorful for medium sized residuals (< 4 and > 2),
  grey/white for small residuals (< 2). More categories or a continuous scale can
  be specified by setting \code{interpolate}. Furthermore, the result of a significance
  test can be visualized by the amount of grey in the colors. If significant, a colorful
  palette is used, if not, the amount of color is reduced.
  See Zeileis, Meyer, and Hornik (2007) and \code{\link[colorspace]{diverge_hcl}} for more details.
  
  The shading \code{shading_max} is applicable in 2-way contingency tables and uses
  a similar strategy as \code{shading_hcl}. But instead of using the cut-offs 2 and 4,
  it employs the critical values for the maximum statistic (by default at 90\% and 99\%).
  Consequently, color in the plot signals a significant result at 90\% or 99\% significance
  level, respectively. The test is carried out by calling \code{\link{coindep_test}}.
  
  The shading \code{shading_Friendly} is very similar to \code{shading_hsv}, but additionally
  codes the sign of the residuals by different line types. See Friendly
  (1994) for more details. \code{shading_Friendly2} and
  \code{shading_sieve} are similar, but use HCL colors.
  
  The shading \code{shading_binary} just visualizes the sign of the residuals by using
  two different colors (default: blue HCL(260, 50, 70) and red HCL(0, 50, 70)).

  \code{shading_Marimekko} is a simple generating function for
  producing, in conjunction with \code{\link{mosaic}}, so-called
  \emph{Marimekko-charts}, which paint the tiles of each columns of a
  mosaic display in the same color to better display departures from
  independence.

  \code{shading_diagonal} generates a color shading for basically square
  matrices (or arrays having the first two dimensons of same length)
  visualizing the diagonal cells, and the off-diagonal cells 1, 2, \dots
  steps removed.
  
  
  
  
  The color implementations employed are \code{\link{hsv}}
  from base R and \code{\link[colorspace]{polarLUV}} from the \pkg{colorspace}
  package, respectively. To transform the HCL coordinates to
  a hexadecimal color string (as returned by \code{hsv}), the function
  \code{\link[colorspace]{hex}} is employed. A convenience wrapper \code{hcl2hex}
  is provided.
}

\references{
  Friendly M. (1994),
  Mosaic Displays for Multi-Way Contingency Tables.
  \emph{Journal of the American Statistical Association},
  \bold{89}, 190--200.

  Meyer D., Zeileis A., and Hornik K. (2006),
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17}(3), 1--48.
  URL http://www.jstatsoft.org/v17/i03/. See also \code{vignette("strucplot", package = "vcd")}.

  Zeileis A., Meyer D., Hornik K. (2007), Residual-Based Shadings for Visualizing
  (Conditional) Independence. \emph{Journal of Computational and Graphical Statistics},
  \bold{16}, 507--525.

  Zeileis A., Hornik K. and Murrell P. (2008),
  Escaping RGBland: Selecting Colors for Statistical Graphics.
  \emph{Computational Statistics & Data Analysis}, Forthcoming.
  Preprint available from \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Zeileis+Hornik+Murrell-2009.pdf}.
}

\value{A shading function which takes only a single argument, interpreted as a
  vector/table of residuals, and returns a \code{"gpar"} object with the
  corresponding vector(s)/table(s) of graphical parameter(s).

}

\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link[colorspace]{hex}},
  \code{\link[colorspace]{polarLUV}},
  \code{\link{hsv}},
  \code{\link{mosaic}},
  \code{\link{assoc}},
  \code{\link{strucplot}},
  \code{\link[colorspace]{diverge_hcl}}
}
\examples{
## load Arthritis data
data("Arthritis")
art <- xtabs(~Treatment + Improved, data = Arthritis)

## plain mosaic display without shading
mosaic(art)

## with shading for independence model
mosaic(art, shade = TRUE)
## which uses the HCL shading
mosaic(art, gp = shading_hcl)
## the residuals are too small to have color,
## hence the cut-offs can be modified
mosaic(art, gp = shading_hcl, gp_args = list(interpolate = c(1, 1.8)))
## the same with the Friendly palette 
## (without significance testing)
mosaic(art, gp = shading_Friendly, gp_args = list(interpolate = c(1, 1.8)))

## assess independence using the maximum statistic
## cut-offs are now critical values for the test statistic
mosaic(art, gp = shading_max)

## association plot with shading as in base R
assoc(art, gp = shading_binary(col = c(1, 2)))

## Marimekko Chart
hec <- margin.table(HairEyeColor, 1:2)
mosaic(hec, gp = shading_Marimekko(hec))
mosaic(HairEyeColor, gp = shading_Marimekko(HairEyeColor))

## Diagonal cells shading
ac <- xtabs(VisualAcuity)
mosaic(ac, gp = shading_diagonal(ac))

}
\keyword{hplot}
