\name{Ord.plot}
\alias{Ord.plot}
\alias{Ord.estimate}
\title{Ord plots}
\description{
Ord plots for diagnosing discrete distributions.
}
\usage{
Ord.plot(obj, legend = TRUE, estimate = TRUE, tol = 0.1, type = NULL,
         ylim = NULL, xlab = "Number of occurrences",
	 ylab = "Frequency ratio", main = "Ord plot", ...)
Ord.estimate(x, type = NULL, tol = 0.1)
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
  \item{ylim}{limits for the y axis.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{main}{a title for the plot.}
  \item{...}{further arguments passed to \code{\link{plot}}.}
  \item{x}{a vector giving intercept and slope for the (fitted)
    line in the Ord plot.}

}
\details{
The Ord plot plots the number of occurrences against a certain
frequency ratio (see Friendly (2000) for details) and
should give a straight line if the data comes from
a poisson, binomial, negative binomial or log-series distribution. The
intercept and slope of this straight line conveys information about
the underlying distribution.

\code{Ord.plot} fits a usual OLS line (black) and a weighted OLS line
(red). From the coefficients of the latter the distribution is
estimated by \code{Ord.estimate} as described in Table 2.10 in Friendly
(2000). To judge whether a coefficient is positive or negative a
tolerance given by \code{tol} is used. If none of the distributions
fits well, no parameters are estimated. Be careful with the conclusions
from \code{Ord.estimate} as it implements just some simple heuristics!
}
\value{a vector giving the intercept and slope of the weighted OLS line.}
\author{
  Achim Zeileis
}
\references{
J. K. Ord (1967), Graphical methods for a class of discrete distributions,
\emph{Journal of the Royal Statistical Society}, \bold{A 130}, 232-238.

Michael Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(1000, size = 1.5, prob = 0.8)
Ord.plot(dummy)

## Real data examples:
data(HorseKicks)
data(Federalist)
data(Butterfly)
data(WomenQueue)
par(mfrow = c(2,2))
Ord.plot(HorseKicks, main = "Death by horse kicks")
Ord.plot(Federalist, main = "Instances of 'may' in Federalist papers")
Ord.plot(Butterfly, main = "Butterfly species collected in Malaya")
Ord.plot(WomenQueue, main = "Women in queues of length 10")
par(mfrow = c(1,1))
}
\keyword{category}
