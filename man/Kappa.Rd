\name{Kappa}
\alias{Kappa}
\alias{print.Kappa}
\alias{summary.Kappa}
\alias{print.summary.Kappa}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cohen's Kappa and weighted Kappa}
\description{
  Computes two agreement rates: Cohen's kappa and weighted kappa, and
  confidence bands.
}
\usage{
Kappa(x, weights = c("Equal-Spacing", "Fleiss-Cohen"), conf.level = 0.95)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a confusion matrix.}
  \item{weights}{either one of the two options or a user-specified
    matrix with same dimensions as \code{x}.}
  \item{conf.level}{level for the confidence intervals.}
}
\details{
  Cohen's kappa is the diagonal sum of the (possibly weighted) relative frequencies,
  corrected for expected values and standardized by its maximum value.
  The equal-spacing weights are defined by \eqn{1 - abs(i - j) / (r -
    1)}, r number of colums/rows, and the Fleiss-Cohen weights by
  \eqn{1 - abs(i - j)^2 / (r - 1)^2}. The latter ones attach greater
  importance to near disagreements.
}
\value{
  An object of class \code{kappa} with three components:
  \item{Kappa}{Kappa statistic, along with Approximate Standard Error
    (ASE) and 95\% confidence bounds.}
  \item{Kappa.Weighted}{idem for the weighted kappa.}
  \item{Weights}{weight matrix used.}
}
\references{
    \item Cohen, Jacob (1960): A coefficient of agreement for nominal
    scales. Educational and Psychological Measurement.
    \item Everitt, B.S. (1968): Moments of statistics kappa and weighted
    kappa. The British Journal of Mathematical and Statistical Psychology.
}
\author{David Meyer\cr
\email{david.meyer@ci.tuwien.ac.at}}
\seealso{\code{\link{agreementplot}}}
\examples{
data(SexualFun)
Kappa(SexualFun)
}
\keyword{category}% at least one, from doc/KEYWORDS
