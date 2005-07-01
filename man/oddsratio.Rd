\name{oddsratio}
\alias{oddsratio}
\alias{print.oddsratio}
\alias{summary.oddsratio}
\alias{confint.oddsratio}
\alias{print.summary.oddsratio}
\alias{plot.oddsratio}
\title{(Log) Odds Ratios}
\description{
  Computes (log) odds ratios and their asymptotic standard errors for
  (possibly) stratified data.
}
\usage{
oddsratio(x, stratum = NULL, log = TRUE)
\method{plot}{oddsratio}(x, conf_level = 0.95, type = "o",
  xlab = "Strata", ylab = NULL, ylim = NULL, whiskers = 0.1, \dots)
}
\arguments{
  \item{x}{a 2 by 2 by \dots table.}
  \item{stratum}{vector of strata dimensions.}
  \item{log}{if \code{FALSE}, ordinary odds ratios are computed.}
  \item{conf_level}{if not \code{NULL} or \code{FALSE},
    \code{conf_level}-\% confidence intervals are plotted for each data
    point.}
  \item{type}{plot type.}
  \item{xlab}{label for the x-axis.}
  \item{ylab}{label for the y-axis.}
  \item{ylim}{y-axis limits.}
  \item{whiskers}{width of the confidence interval whiskers.}
  \item{\dots}{other graphics parameters (see \code{\link{par}}).}
}
\value{
  An object of class \code{"logoddsratio"}, which is simply a vector of
  (log) odds ratios with dimensionality depending on \code{stratum},
  along with the following attributes:
  \item{ASE}{a numeric vector with the asymptotic standard errors.}
  \item{log}{logical indicating whether log odds ratios or common odds
    ratios are computed.}
}
\note{
  The \code{summary} method prints the standard errors and---for log
  odds ratios---also computes and prints asymptotic \eqn{z} tests
  (standardized log odds ratios) and the corresponding \eqn{p} values.

  There is a \code{confint} method for computing confidence intervals
  for the (log) odds ratios.

  The \code{plot} method plots (log) odds ratios, computed by \code{oddsratio}
  for \eqn{2 \times 2 \times k}{2 x 2 x k} tables, along with confidence
  intervals.
}

\references{
  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\seealso{
  \code{\link{confint}}
}
\examples{
## load Coal Miners data
data(CoalMiners)

## compute log odds ratios
lor <- oddsratio(CoalMiners)
lor

## summary with z tests
summary(lor)

## confidence intervals
confint(lor)

## visualization
plot(lor,
     xlab = "Age Group",
     main = "Breathelessness and Wheeze in Coal Miners")

## add quadratic model
g <- seq(25, 60, by = 5)
m <- lm(lor ~ g + I(g^2))
lines(fitted(m), col = "red")
}
\keyword{category}
