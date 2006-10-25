\name{legends}
\alias{legends}
\alias{legend_resbased}
\alias{legend_fixed}
\title{Legend Functions for Strucplots}
\description{
  These functions generate legend functions for residual-based shadings.
}
\usage{
legend_resbased(fontsize = 12,
  x = unit(1, "lines"), y = unit(0.1,"npc"),
  height = unit(0.8, "npc"),
  width = unit(0.7, "lines"),
  digits = 3, check_overlap = TRUE, text = NULL,
  steps = 200, ticks = 10, pvalue = TRUE)
legend_fixed(fontsize = 12, x = unit(1, "lines"), y = NULL,
  height = NULL, width = unit(1.5, "lines"), steps = 200,
  digits = 1, space = 0.05, text = NULL)
}
\arguments{
  \item{fontsize}{fontsize of title and p-value text.}
  \item{x, y}{objects of class \code{"unit"} indicating the coordinates
    of the title. For \code{legend_fixed}, the default for \code{y} is
    computed as to leave enough space for the specified \code{text}.}
  \item{height, width}{object of class \code{"unit"} indicating the
    height/width of the legend. For \code{legend_fixed}, the default for
    \code{y} is computed as to align upper margins of legend and actual plot.} 
  \item{digits}{number of digits for the scale labels.}
  \item{check_overlap}{logical indicating whether overlap of scale
    labels should be inhibited or not.}
  \item{space}{For \code{legend_fixed} only: proportion of space
    between the tiles.}
  \item{text}{character string indicating the title of the legend.}
  \item{steps}{granularity of the color gradient.}
  \item{ticks}{number of scale ticks.}
  \item{pvalue}{logical indicating whether the \eqn{p}-value should be
    visualized under the scale or not.}
}
\value{
  A function with arguments:
  \item{residuals}{residuals from the fitted independence model to be
    visualized.}
  \item{shading}{shading function computing colors from residuals (see
    details).}
  \item{autotext}{character vector indicating the title to be used when
    no \code{text} argument is specified. Allows strucplot to generate
    sensible defaults depending on the residuals type.}
}
\details{
  These functions generate legend functions for residual-based shadings,
  visualizing deviations from expected values of an hypothesized
  independence model.  Therefore, the legend uses a supplied shading
  function to visualize the color gradient for the residuals range.
  \code{legend_fixed} is inspired by the legend used in
  \code{\link[graphics]{mosaicplot}}. For
  more details on the shading functions and their return values, see
  \code{\link{shadings}}.
}
\references{
  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.

  Meyer, D., Zeileis, A., Hornik, K. (2003), Visualizing independence using
  extended association plots. \emph{Proceedings of the 3rd International
  Workshop on Distributed Statistical Computing}, K. Hornik, F. Leisch,
  A. Zeileis (eds.), ISSN 1609-395X.
  \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
}

\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{structable}},
  \code{\link{shadings}}
}

\examples{
data("Titanic")

mosaic(Titanic, shade = TRUE, legend = legend_resbased)
mosaic(Titanic, shade = TRUE, legend = legend_fixed, gp = shading_Friendly)
}
\keyword{hplot}

