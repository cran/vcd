\name{agreementplot.default}
\alias{agreementplot}
\alias{agreementplot.default}
\alias{agreementplot.formula}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Bangdiwala's Observer Agreement Chart}
\description{
  Representation of a k x k confusion matrix, where the
  observed and expected diagonal elements are represented by superposed
  black and white rectangles, respectively. The function also computes
  a statistic measuring the strength of agreement (relation of respective
  area sums).
}
\usage{
\method{agreementplot}{default}(x, reverse.y = TRUE, main = "Agreement Chart",
              weights = c(1, 1 - 1/(ncol(x) - 1)^2), cex.main = 2, cex.lab = 1.5,
              xlab = names(dimnames(x))[2], ylab = names(dimnames(x))[1], ...)
\method{agreementplot}{formula}(formula, data = NULL, ..., subset)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a confusion matrix, i.e. a table with equal-sized dimensions.}
  \item{reverse.y}{if \code{TRUE}, the y axis is reversed (i.e., the
    rectangles' positions correspond to the contingency table).}
  \item{main}{user-specified main title.}
  \item{weights}{vector of weights for successive larger observed areas, used
    in the agreement strength statistic, and also for the shading. The
    first element should be 1.}
  \item{cex.main}{font size of title.}
  \item{cex.lab}{font size of labels.}
  \item{xlab, ylab}{labels of x- and y-axis.}
  \item{formula}{a formula, such as \code{y ~ x}. For details, see \code{xtabs}.}
  \item{data}{a data.frame (or list), or a contingency table from which the
          variables in `formula' should be taken.}
  \item{subset}{an optional vector specifying a subset of observations to be
    used for plotting.}
  \item{\dots}{further graphics parameters (see \code{par}).}
}
\details{
  Weights can be specified to allow for partial agreement, taking into
  account contributions from off-diagonal cells. A weight vector of
  length 1 means strict agreement only, each additional element
  increases the maximum number of disagreement steps.
}
\value{
  Invisibly returned, a list with components
  \item{Bangdiwala}{the unweighted agreement strength statistic}
  \item{Bangdiwala.Weighted}{the weighted statistic}
  \item{weights}{the weigtht vector used.}
}
\references{
Michael Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}

\examples{
data(SexualFun)
agreementplot(t(SexualFun))

data(MSPatients)
## Enlarge plot manually or use sth. like: X11(width = 12)!
par(mfrow = c(1,2))
agreementplot(t(MSPatients[,,1]), main = "Winnipeg Patients")
agreementplot(t(MSPatients[,,2]), main = "New Orleans Patients")
par(mfrow = c(1,1))
}
\keyword{category}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
