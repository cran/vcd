\name{assocstats}
\alias{assocstats}
\alias{summary.assocstats}
\alias{print.assocstats}
\alias{print.summary.assocstats}
\title{Association Statistics}
\description{
  Computes the Pearson chi-Squared test, the Likelihood Ratio chi-Squared test,
  the phi coefficient, the contingency coefficient and Cramer's V for
  possibly stratified contingency tables.
}
\usage{
assocstats(x)
}
\arguments{
  \item{x}{a contingency table, with possibly more than 2 dimensions. In
  this case, all dimensions except the first two ones are considered as strata.}
}
\value{
  In case of a 2-dimensional table, a list with components:
  \item{chisq_tests}{a \eqn{2 \times 3}{2 x 3} table with the
    chi-squared statistics.}
  \item{phi}{The \emph{absolute value} of the phi coefficient (only
    defined for \eqn{2 \times 2}{2 x 2} tables).}
  \item{cont}{The contingency coefficient.}
  \item{cramer}{Cramer's V.}
  In case of higher-dimensional tables, a list of the above mentioned
  structure, each list component representing one stratum defined by the
  combinations of all levels of the stratum dimensions.
}
\references{
  Michael Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.

  Fleiss, J. L. (1981).
  \emph{Statistical methods for rates and proportions}
  (2nd ed). New York: Wiley 
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data("Arthritis")
tab <- xtabs(~Improved + Treatment, data = Arthritis)
summary(assocstats(tab))

assocstats(UCBAdmissions)
}
\keyword{category}
