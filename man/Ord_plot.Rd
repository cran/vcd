\name{Ord_plot}
\alias{Ord_plot}
\alias{Ord_estimate}
\title{Ord Plots}
\description{
  Ord plots for diagnosing discrete distributions.
}
\usage{
Ord_plot(obj, legend = TRUE, estimate = TRUE, tol = 0.1, type = NULL,
  xlim = NULL, ylim = NULL, xlab = "Number of occurrences",
  ylab = "Frequency ratio", main = "Ord plot", gp = gpar(cex = 0.5),
  name = "Ord_plot", newpage = TRUE, pop = TRUE, ...)
Ord_estimate(x, type = NULL, tol = 0.1)
}
\arguments{
  \item{obj}{either a vector of counts, a 1-way table of frequencies
    of counts or a data frame or matrix with frequencies
    in the first column and the corresponding counts in the
    second column.}
  \item{legend}{logical. Should a legend be plotted?.}
  \item{estimate}{logical. Should the distribution and its parameters
    be estimated from the data? See details.}
  \item{tol}{tolerance for estimating the distribution. See details.}
  \item{type}{a character string indicating the distribution, must be
    one of \code{"poisson"}, \code{"binomial"}, \code{"nbinomial"}
    or \code{"log-series"} or \code{NULL}. In the latter case the
    distribution is estimated from the data. See details.}
  \item{xlim}{limits for the x axis.}
  \item{ylim}{limits for the y axis.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{main}{a title for the plot.}
  \item{gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the points.}  
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
  \item{...}{further arguments passed to \code{\link{grid.points}}.}
  \item{x}{a vector giving intercept and slope for the (fitted)
    line in the Ord plot.}

}
\details{
  The Ord plot plots the number of occurrences against a certain
  frequency ratio (see Friendly (2000) for details) and should give a
  straight line if the data comes from a poisson, binomial, negative
  binomial or log-series distribution.  The intercept and slope of this
  straight line conveys information about the underlying distribution.

  \code{Ord_plot} fits a usual OLS line (black) and a weighted OLS line
  (red).  From the coefficients of the latter the distribution is
  estimated by \code{Ord_estimate} as described in Table 2.10 in
  Friendly (2000).  To judge whether a coefficient is positive or
  negative a tolerance given by \code{tol} is used.  If none of the
  distributions fits well, no parameters are estimated.  Be careful with
  the conclusions from \code{Ord_estimate} as it implements just some
  simple heuristics!
}
\value{
  A vector giving the intercept and slope of the weighted OLS line.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\references{
  J. K. Ord (1967),
  Graphical methods for a class of discrete distributions,
  \emph{Journal of the Royal Statistical Society}, \bold{A 130},
  232--238.

  Michael Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(1000, size = 1.5, prob = 0.8)
Ord_plot(dummy)

## Real data examples:
data(HorseKicks)
data(Federalist)
data(Butterfly)
data(WomenQueue)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
Ord_plot(HorseKicks, main = "Death by horse kicks", newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
Ord_plot(Federalist, main = "Instances of 'may' in Federalist papers", newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
Ord_plot(Butterfly, main = "Butterfly species collected in Malaya", newpage = FALSE)
popViewport()

pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
Ord_plot(WomenQueue, main = "Women in queues of length 10", newpage = FALSE)
popViewport(2)
}
\keyword{category}
