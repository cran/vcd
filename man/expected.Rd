\name{expected}
\alias{expected}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Expected Frequencies Table}
\description{
Computes table of expected frequencies from a n-way table.
}
\usage{
expected(x, frequency = c("absolute", "relative"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a table.}
  \item{frequency}{indicates, if absolute or relative frequencies should
  be computed.}
}
\value{
  A table with either absolute or relative frequencies.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\examples{
data(MSPatients)
expected(MSPatients)
expected(MSPatients, frequency = "relative")
}
\keyword{category}% at least one, from doc/KEYWORDS
\keyword{array}% __ONLY ONE__ keyword per line
