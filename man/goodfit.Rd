\name{goodfit}
\alias{goodfit}
\alias{summary.goodfit}
\alias{plot.goodfit}
\alias{predict.goodfit}
\alias{fitted.goodfit}
\alias{print.goodfit}
\title{Goodness-of-fit Tests for Discrete Data}
\description{
  Fits a discrete (count data) distribution for goodness-of-fit tests.
}
\usage{
goodfit(x, type = c("poisson", "binomial", "nbinomial"),
  method = c("ML", "MinChisq"), par = NULL)
\method{predict}{goodfit}(object, newcount = NULL, type = c("response", "prob"), \dots)
}
\arguments{
  \item{x}{either a vector of counts, a 1-way table of frequencies of
    counts or a data frame or matrix with frequencies in the first
    column and the corresponding counts in the second column.}
  \item{type}{a character string indicating which distribution should be
    fit (for \code{goodfit}) or  indicating the type of prediction
    (fitted response or probabilities in \code{predict}) respectively.}
  \item{method}{a character string indicating whether the distribution
    should be fit via ML (Maximum Likelihood) or Minimum Chi-squared.}
  \item{par}{a named list giving the distribution parameters (named as
    in the corresponding density function), if set to \code{NULL}, the
    default, the parameters are estimated.  If the parameter \code{size}
    is not specified if \code{type} is \code{"binomial"} it is taken to
    be the maximum count. If \code{type} is \code{"nbinomial"}, then
    parameter \code{size} can be specified to fix it so that only the
    parameter \code{prob} will be estimated (see the examples below).}
  \item{object}{an object of class \code{"goodfit"}.}
  \item{newcount}{a vector of counts.  By default the counts stored in
    \code{object} are used, i.e., the fitted values are computed.  These
    can also be extracted by \code{fitted(object)}.}
  \item{\dots}{\emph{currently not used}.}
}
\details{
  \code{goodfit} essentially computes the fitted values of a discrete
  distribution (either poisson, binomial or negative binomial) to the
  count data given in \code{x}.  If the parameters are not specified
  they are estimated either by ML or Minimum Chi-squared. 

  \code{par} should be a named list specifying the parameters \code{lambda}
  for \code{"poisson"} and \code{prob} and \code{size} for
  \code{"binomial"} or \code{"nbinomial"}, respectively.
  If for \code{"binomial"} \code{size} is not specified it is not
  estimated but taken as the maximum count.

  The corresponding Pearson Chi-squared or likelihood ratio statistic
  respectively is computed and given with their \eqn{p} values by the
  \code{summary} method.  The \code{plot} method produces a
  \code{\link{rootogram}} of the observed and fitted values.
}
\value{
  A list of class \code{"goodfit"} with elements:
  \item{observed}{observed frequencies.}
  \item{count}{corresponding counts.}
  \item{fitted}{expected frequencies (fitted by ML).}
  \item{type}{a character string indicating the distribution fitted.}
  \item{method}{a character string indicating the fitting method (can
    be either \code{"ML"}, \code{"MinChisq"} or \code{"fixed"} if the
    parameters were specified).}
  \item{df}{degrees of freedom.}
  \item{par}{a named list of the (estimated) distribution parameters.}
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\references{
  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(200, size = 1.5, prob = 0.8)
gf <- goodfit(dummy, type = "nbinomial", method = "MinChisq")
summary(gf)
plot(gf)

dummy <- rbinom(100, size = 6, prob = 0.5)
gf1 <- goodfit(dummy, type = "binomial", par = list(size = 6))
gf2 <- goodfit(dummy, type = "binomial", par = list(prob = 0.6, size = 6))
summary(gf1)
plot(gf1)
summary(gf2)
plot(gf2)

## Real data examples:
data("HorseKicks")
HK.fit <- goodfit(HorseKicks)
summary(HK.fit)
plot(HK.fit)

data("Federalist")
## try geometric and full negative binomial distribution
F.fit <- goodfit(Federalist, type = "nbinomial", par = list(size = 1))
F.fit2 <- goodfit(Federalist, type = "nbinomial")
summary(F.fit)
summary(F.fit2)
plot(F.fit)
plot(F.fit2)
}
\keyword{category}
