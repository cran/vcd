\name{oddsratio}
\alias{oddsratio}
\alias{print.oddsratio}
\alias{summary.oddsratio}
\alias{print.summary.oddsratio}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Odds Ratios}
\description{
  Computes (log) odds ratios and their asymtotic standard errors
  for (possibly) stratified data.
}
\usage{
oddsratio(x, stratum = NULL, log = TRUE, conf.level = 0.95)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a 2 x 2 x ... table.}
 \item{stratum}{vector of strata dimensions.}
 \item{log}{if \code{FALSE}, ordinary odds ratios are computed}
 \item{conf.level}{confidence level used for confidence intervals.}
}
\value{
  An object of class \code{logoddratio}, which is simply a
  vector of (log) odds ratios with dimensionality depending on \code{stratum}, along with
  the following attributes:
  \item{ASE}{a vector with the asymptotic standard errors}
  \item{lwr, upr}{a vector with lower/upper confidence bounds}
  \item{Z}{test statistic for significance tests (only computed for log
  odds ratios: simply the standardized log odds ratio, which follows
an asymptotic normal distribution.}
  \item{P}{p-value for this test.}
  \item{log}{logical value indicating if log odds ratios or common odds
    ratios are computed.}
}
\note{The \code{summary} method prints the standard errors.}

\references{
M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}

\seealso{\code{\link{plot.oddsratio}}}

\examples{
data(CoalMiners)
summary(oddsratio(CoalMiners))
}
\keyword{category}% __ONLY ONE__ keyword per line
