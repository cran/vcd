\name{independence_table}
\alias{independence_table}
\title{Independence Table}
\description{
  Computes table of expected frequencies (under the null hypotheses of
  independence) from an \eqn{n}-way table.
}
\usage{
independence_table(x, frequency = c("absolute", "relative"))
}
\arguments{
  \item{x}{a table.}
  \item{frequency}{indicates whether absolute or relative frequencies
    should be computed.}
}
\value{
  A table with either absolute or relative frequencies.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data(MSPatients)
independence_table(MSPatients)
independence_table(MSPatients, frequency = "relative")
}
\keyword{category}
\keyword{array}
