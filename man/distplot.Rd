\name{distplot}
\alias{distplot}
\title{Diagnostic distribution plots}
\description{
Diagnostic distribution plots: poissonness, binomialness and
negative bionomialness plots.
}
\usage{
distplot(obj, type = c("poisson", "binomial", "nbinomial"),
         size = NULL, lambda = NULL, legend = TRUE, ylim = NULL,
         line.col = 2, conf.int = TRUE, conf.level = 0.95, main = NULL,
	 xlab = "Number of occurrences", ylab = "Distribution metameter", ...)
}
\arguments{
  \item{obj}{either a vector of counts, a 1-way table of frequencies
    of counts or a data frame or matrix with frequencies
    in the first column and the corresponding counts in the
    second column.}
  \item{type}{a character string indicating the distribution.}
  \item{size}{the size argument for the binomial distribution, only
    required if \code{type} is \code{"binomial"}. If set to \code{NULL},
    \code{size} is taken to be the maximum count.}
  \item{lambda}{parameter of the poisson distribution. If type is
    \code{"poisson"} and \code{lambda} is specified a leveled
    poissonness plot is produced.}
  \item{legend}{logical. Should a legend be plotted?.}
  \item{ylim}{limits for the y axis.}
  \item{line.col}{color for fitted line.}
  \item{conf.int}{logical. Should confidence intervals be plotted?.}
  \item{conf.level}{confidence level for confidence intervals.}
  \item{main}{a title for the plot.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{...}{further arguments passed to \code{\link{plot}}.}
}
\details{
\code{distplot} plots the number of occurrences (counts) against the
distribution metameter of the specified distribution. If
the distribution fits the data, the plot should show a
straight line. See Friendly (2000) for details.
}
\author{
  Achim Zeileis
}
\references{
D. C. Hoaglin (1980), A poissonness plot, \emph{The American
Statistican}, \bold{34}, 146-149.

D. C. Hoaglin & J. W. Tukey (1985), Checking the shape of discrete
distributions. In D. C. Hoaglin, F. Mosteller, J. W. Tukey (eds.),
\emph{Exploring Data Tables, Trends and Shapes}, chapter 9. John
Wiley & Sons, New York.

Michael Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(1000, size = 1.5, prob = 0.8)
distplot(dummy, type = "nbinomial")

## Real data examples:
data(HorseKicks)
data(Federalist)
data(Saxony)
distplot(HorseKicks, type = "poisson")
distplot(HorseKicks, type = "poisson", lambda = 0.61)
distplot(Federalist, type = "poisson")
distplot(Federalist, type = "nbinomial")
distplot(Saxony, type = "binomial", size = 12)
}
\keyword{category}
