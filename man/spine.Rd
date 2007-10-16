\name{spine}
\alias{spine}
\alias{spine.default}
\alias{spine.formula}
\title{Spine Plots and Spinograms}
\description{
  Spine plots are a special cases of mosaic plots, and can be seen as
  a generalization of stacked (or highlighted) bar plots. Analogously,
  spinograms are an extension of histograms.
}
\usage{
spine(x, \dots)
\method{spine}{default}(x, y = NULL,
  breaks = NULL, ylab_tol = 0.05, off = NULL,
  main = "", xlab = NULL, ylab = NULL, ylim = c(0, 1), margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "spineplot", newpage = TRUE, pop = TRUE,
  \dots)
\method{spine}{formula}(formula, data = list(),
  breaks = NULL, ylab_tol = 0.05, off = NULL,
  main = "", xlab = NULL, ylab = NULL, ylim = c(0, 1), margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "spineplot", newpage = TRUE, pop = TRUE,
  \dots)
}
\arguments{
  \item{x}{an object, the default method expects either a single variable
    (interpreted to be the explanatory variable) or a 2-way table. See details.}
  \item{y}{a \code{"factor"} interpreted to be the dependent variable}
  \item{formula}{a \code{"formula"} of type \code{y ~ x} with a single dependent \code{"factor"}
    and a single explanatory variable.}    
  \item{data}{an optional data frame.}
  \item{breaks}{if the explantory variable is numeric, this controls how
    it is discretized. \code{breaks} is passed to \code{\link{hist}} and can
    be a list of arguments.}
  \item{ylab_tol}{convenience tolerance parameter for y-axis annotation.
    If the distance between two labels drops under this threshold, they are
    plotted equidistantly.}
  \item{off}{vertical offset between the bars (in per cent). It is fixed to
    \code{0} for spinograms and defaults to \code{2} for spine plots.}
  \item{main, xlab, ylab}{character strings for annotation}
  \item{ylim}{limits for the y axis}
  \item{margins}{margins when calling \code{\link{plotViewport}}}
  \item{gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the rectangles. It should specify in particular a vector of
    \code{fill} colors of the same length as \code{levels(y)}. The default is
    to call \code{\link{gray.colors}}.}  
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
  \item{\dots}{additional arguments passed to \code{\link{plotViewport}}.}
}
\details{
  \code{spine} creates either a spinogram or a spine plot. It can be called
  via \code{spine(x, y)} or \code{spine(y ~ x)} where \code{y} is interpreted
  to be the dependent variable (and has to be categorical) and \code{x}
  the explanatory variable. \code{x} can be either categorical (then a spine
  plot is created) or numerical (then a spinogram is plotted).
  Additionally, \code{spine} can also be called with only a single argument
  which then has to be a 2-way table, interpreted to correspond to \code{table(x, y)}.
  
  Spine plots are a generalization of stacked bar plots where not the heights
  but the widths of the bars corresponds to the relative frequencies of \code{x}.
  The heights of the bars then correspond to the conditional relative frequencies
  of \code{y} in every \code{x} group. This is a special case of a mosaic plot
  with specific spacing and shading.
  
  Analogously, spinograms extend stacked histograms. As for the histogram,
  \code{x} is first discretized (using \code{\link{hist}}) and then for the
  discretized data a spine plot is created.
}
\value{
  The table visualized is returned invisibly.
}
\seealso{
  \code{\link{cd_plot}}, \code{\link{mosaic}}, \code{\link{hist}}
}
\references{
  Hummel, J. (1996), Linked bar charts: Analysing categorical data graphically.
  \emph{Computational Statistics}, \bold{11}, 23--33.

  Hofmann, H., Theus, M. (2005), \emph{Interactive graphics for visualizing
  conditional distributions}, Unpublished Manuscript.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
## Arthritis data (dependence on a categorical variable)
data("Arthritis")
(spine(Improved ~ Treatment, data = Arthritis))

## Arthritis data (dependence on a numerical variable)
(spine(Improved ~ Age, data = Arthritis, breaks = 5))
(spine(Improved ~ Age, data = Arthritis, breaks = quantile(Arthritis$Age)))
(spine(Improved ~ Age, data = Arthritis, breaks = "Scott"))

## Space shuttle data (dependence on a numerical variable)
data("SpaceShuttle")
(spine(Fail ~ Temperature, data = SpaceShuttle, breaks = 3))
}
\keyword{hplot}
