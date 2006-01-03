\name{assocstats}
\alias{assocstats}
\alias{summary.assocstats}
\alias{print.assocstats}
\alias{print.summary.assocstats}
\title{Association Statistics}
\description{
  Computes the Pearson chi-Squared test, the Likelihood Ratio chi-Squared test,
  the phi coefficient, the contingency coefficient and Cramer's V.
}
\usage{
assocstats(x)
}
\arguments{
  \item{x}{an \eqn{r \times c}{r x c} table.}
}
\value{
  A list with components:
  \item{chisq.tests}{a \eqn{2 \times 3}{2 x 3} table with the
    chi-squared statistics.}
  \item{phi}{The phi coefficient.}
  \item{cont}{The contingency coefficient.}
  \item{cramer}{Cramer's V.}
}
\references{
  Michael Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data("Arthritis")
tab <- xtabs(~Improved + Treatment, data = Arthritis)
summary(assocstats(tab))
}
\keyword{category}
