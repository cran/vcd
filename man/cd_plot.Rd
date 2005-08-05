\name{cd_plot}
\alias{cd_plot}
\alias{cd_plot.default}
\alias{cd_plot.formula}
\title{Conditional Density Plots}
\description{
  Computes and plots conditional densities describing how the
  distribution of a categorical variable \code{y} changes over a
  numerical variable \code{x}.
}
\usage{
cd_plot(x, \dots)
\method{cd_plot}{default}(x, y,
  plot = TRUE, ylab_tol = 0.05,
  bw = "nrd0", n = 512, from = NULL, to = NULL,
  main = "", xlab = NULL, ylab = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "cd_plot", newpage = TRUE, pop = TRUE,
  \dots)
\method{cd_plot}{formula}(formula, data = list(),
  plot = TRUE, ylab_tol = 0.05,
  bw = "nrd0", n = 512, from = NULL, to = NULL,
  main = "", xlab = NULL, ylab = NULL, margins = c(5.1, 4.1, 4.1, 3.1),
  gp = gpar(), name = "cd_plot", newpage = TRUE, pop = TRUE,
  \dots)
}
\arguments{
  \item{x}{an object, the default method expects either a single numerical variable.}
  \item{y}{a \code{"factor"} interpreted to be the dependent variable}
  \item{formula}{a \code{"formula"} of type \code{y ~ x} with a single dependent \code{"factor"}
    and a single numerical explanatory variable.}    
  \item{data}{an optional data frame.}
  \item{plot}{logical. Should the computed conditional densities be plotted?}
  \item{ylab_tol}{convenience tolerance parameter for y-axis annotation.
    If the distance between two labels drops under this threshold, they are
    plotted equidistantly.}
  \item{bw, n, from, to, \dots}{arguments passed to \code{\link{density}}}
  \item{main, xlab, ylab}{character strings for annotation}
  \item{margins}{margins when calling \code{\link{plotViewport}}}
  \item{gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the rectangles. It should specify in particular a vector of
    \code{fill} colors of the same length as \code{levels(y)}. The default is
    to call \code{\link{gray.colors}}.}  
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
}
\details{
  \code{cd_plot} computes the conditional densities of \code{x} given
  the levels of \code{y} weighted by the marginal distribution of \code{y}.
  The densities are derived cumulatively over the levels of \code{y}.

  This visualization technique is similar to spinograms (see \code{\link{spine}})
  but they do not discretize the explanatory variable, but rather use a smoothing
  approach. Furthermore, the original x axis and not a distorted x axis (as for
  spinograms) is used. This typically results in conditional densities that
  are based on very few observations in the margins: hence, the estimates are less
  reliable there.
}
\value{
  The conditional density functions (cumulative over the levels of \code{y})
  are returned invisibly.
}
\seealso{
  \code{\link{spine}}, \code{\link{density}}
}
\references{
  Hofmann, H., Theus, M. (2005), \emph{Interactive graphics for visualizing
  conditional distributions}, Unpublished Manuscript.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
## Arthritis data
data("Arthritis")
cd_plot(Improved ~ Age, data = Arthritis)
cd_plot(Improved ~ Age, data = Arthritis, bw = 3)
cd_plot(Improved ~ Age, data = Arthritis, bw = "SJ")
## compare with spinogram
spine(Improved ~ Age, data = Arthritis, breaks = 3)

## Space shuttle data
data("SpaceShuttle")
cd_plot(Fail ~ Temperature, data = SpaceShuttle, bw = 2)

## scatter plot with conditional density
cdens <- cd_plot(Fail ~ Temperature, data = SpaceShuttle, bw = 2, plot = FALSE)
plot(I(-1 * (as.numeric(Fail) - 2)) ~ jitter(Temperature, factor = 2), data = SpaceShuttle,
  xlab = "Temperature", ylab = "Failure")
lines(53:81, cdens[[1]](53:81), col = 2)
}
\keyword{hplot}
