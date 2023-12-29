\name{Kappa}
\alias{Kappa}
\alias{print.Kappa}
\alias{confint.Kappa}
\alias{summary.Kappa}
\alias{print.summary.Kappa}
\title{Cohen's Kappa and Weighted Kappa}
\description{
  Computes two agreement rates: Cohen's kappa and weighted kappa, and
  confidence bands.
}
\usage{
Kappa(x, weights = c("Equal-Spacing", "Fleiss-Cohen"))
\S3method{print}{Kappa}(x, digits=max(getOption("digits") - 3, 3),
CI=FALSE, level=0.95, ...)
\S3method{confint}{Kappa}(object, parm, level = 0.95, ...)
\S3method{summary}{Kappa}(object, ...)
\S3method{print}{summary.Kappa}(x, ...)
}
\arguments{
  \item{x}{For \code{Kappa}: a confusion matrix. For the print methods:
    object of class \code{"Kappa"} or \code{"summary.Kappa"}}
  \item{weights}{either one of the character strings given in the
    default value, or a user-specified matrix with same dimensions as
    \code{x}.}
  \item{digits}{minimal number of significant digits.}
  \item{CI}{logical; shall confidence limits be added to the output?}
  \item{level}{confidence level between 0 and 1 used for the confidence
    interval.}
  \item{object}{object of class \code{"Kappa"}.}
  \item{parm}{Currently, ignored.}
  \item{\dots}{Further arguments passed to the default print method.}
}
\details{
  Cohen's kappa is the diagonal sum of the (possibly weighted) relative
  frequencies, corrected for expected values and standardized by its
  maximum value.
  The equal-spacing weights are defined by \eqn{1 - |i - j| / (r -
    1)}{1 - abs(i - j) / (r - 1)}, \eqn{r} number of columns/rows, and
  the Fleiss-Cohen weights by \eqn{1 - |i - j|^2 / (r - 1)^2}{1 -
    abs(i - j)^2 / (r - 1)^2}.
  The latter one attaches greater importance to near disagreements.
    
}
\value{
  An object of class \code{"Kappa"} with three components:
  \item{Unweighted}{numeric vector of length 2 with the kappa statistic
    (\code{value} component), along with Approximate Standard Error
    (\code{ASE} component)}
  \item{Weighted}{idem for the weighted kappa.}
  \item{Weights}{numeric matrix with weights used.}
}
\note{
  The \code{summary} method also prints the weights.

  There is a \code{confint} method for computing approximate confidence
  intervals.
}
\references{
  Cohen, J. (1960),
  A coefficient of agreement for nominal scales.
  \emph{Educational and Psychological Measurement},
  \bold{20}, 37--46.

  Everitt, B.S. (1968),
  Moments of statistics kappa and weighted kappa.
  \emph{The British Journal of Mathematical and Statistical Psychology},
  \bold{21}, 97--103.

  Fleiss, J.L., Cohen, J., and Everitt, B.S. (1969),
  Large sample standard errors of kappa and weighted kappa.
  \emph{Psychological Bulletin},
  \bold{72}, 332--327.

}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{agreementplot}},
  \code{\link{confint}}
}
\examples{
data("SexualFun")
K <- Kappa(SexualFun)
K
confint(K)
summary(K)
print(K, CI = TRUE)

}
\keyword{category}
