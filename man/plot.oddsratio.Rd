\name{plot.oddsratio}
\alias{plot.oddsratio}
\title{(Log) Odds Ratios Plots}
\description{
Plots (log) odds ratios, computed by \code{oddsratio} for 2 x 2 x
k-tables, along with confidence intervals.
}
\usage{
plot.oddsratio(x, confidence = TRUE,
                  type ="o", ylab = "Log Odds Ratio", xlab = "Strata",
                  whiskers = 0.1, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{logoddsratio}, computed for a 2 x 2
    x k-table.}
  \item{confidence}{if \code{TRUE}, confidence intervals are plotted for
  each data point.}
  \item{type}{plot type.}
  \item{ylab}{label for the y-axis.}
  \item{xlab}{label for the x-axis.}
  \item{whiskers}{Width of the confidence interval whiskers.}
  \item{\dots}{Other graphics parameters (see \code{par}).}
}
\references{
M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}

\seealso{\code{\link{oddsratio}}}

\examples{
data(CoalMiners)
l <- oddsratio(CoalMiners)
g <- seq(25, 60, by = 5)
plot(l,
     xlab = "Age Group",
     main = "Breathelessness and Wheeze in Coal Miners")
m <- lm(l ~ g + I(g^2))
lines(fitted(m), col = "red")
}

\keyword{hplot}% at least one, from doc/KEYWORDS
