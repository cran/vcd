\name{distplot}
\alias{distplot}
\title{Diagnostic Distribution Plots}
\description{
  Diagnostic distribution plots: poissonness, binomialness and
  negative binomialness plots.
}
\usage{
distplot(x, type = c("poisson", "binomial", "nbinomial"),
  size = NULL, lambda = NULL, legend = TRUE, xlim = NULL, ylim = NULL,
  conf_int = TRUE, conf_level = 0.95, main = NULL,
  xlab = "Number of occurrences", ylab = "Distribution metameter",
  gp = gpar(cex = 0.8), lwd=2, name = "distplot", newpage = TRUE, pop = TRUE, \dots)
}
\arguments{
  \item{x}{either a vector of counts, a 1-way table of frequencies of
    counts or a data frame or matrix with frequencies in the first
    column and the corresponding counts in the second column.}
  \item{type}{a character string indicating the distribution.}
  \item{size}{the size argument for the binomial and negative binomial
    distribution.
    If set to \code{NULL} and \code{type} is \code{"binomial"}, then
    \code{size} is taken to be the maximum count. 
    If set to \code{NULL} and \code{type} is \code{"nbinomial"}, then
    \code{size} is estimated from the data.}
  \item{lambda}{parameter of the poisson distribution.
    If type is \code{"poisson"} and \code{lambda} is specified a leveled
    poissonness plot is produced.}
  \item{legend}{logical.  Should a legend be plotted?}
  \item{xlim}{limits for the x axis.}
  \item{ylim}{limits for the y axis.}
  \item{conf_int}{logical.  Should confidence intervals be plotted?}
  \item{conf_level}{confidence level for confidence intervals.}
  \item{main}{a title for the plot.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the points.}  
  \item{lwd}{line width for the fitted line}
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
  \item{...}{further arguments passed to \code{\link{grid.points}}.}
}
\details{
  \code{distplot} plots the number of occurrences (counts) against the
  distribution metameter of the specified distribution.  If the
  distribution fits the data, the plot should show a straight line.
  See Friendly (2000) for details.
  
  In these plots, the open points show the observed count metameters;
  the filled points show the confidence interval centers, and the
  dashed lines show the \code{conf_level} confidence intervals for
  each point.
}

\value{
Returns invisibly a data frame containing the counts (\code{Counts}), 
frequencies (\code{Freq}) and other details of the computations used
to construct the plot.
}

\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\references{
  D. C. Hoaglin (1980),
  A poissonness plot,
  \emph{The American Statistican}, \bold{34}, 146--149.

  D. C. Hoaglin & J. W. Tukey (1985),
  Checking the shape of discrete distributions.
  In D. C. Hoaglin, F. Mosteller, J. W. Tukey (eds.),
  \emph{Exploring Data Tables, Trends and Shapes}, chapter 9.
  John Wiley & Sons, New York.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(1000, size = 1.5, prob = 0.8)
distplot(dummy, type = "nbinomial")

## Real data examples:
data("HorseKicks")
data("Federalist")
data("Saxony")
distplot(HorseKicks, type = "poisson")
distplot(HorseKicks, type = "poisson", lambda = 0.61)
distplot(Federalist, type = "poisson")
distplot(Federalist, type = "nbinomial", size = 1)
distplot(Federalist, type = "nbinomial")
distplot(Saxony, type = "binomial", size = 12)
}
\keyword{category}
