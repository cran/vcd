\name{mar.table}
\alias{mar.table}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Table with marginal sums.}
\description{
  Adds row and column sums to a two-way table.
}
\usage{
mar.table(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a two-way table.}
}
\value{
  A table with added total column and total row.
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\examples{
data(SexualFun)
mar.table(SexualFun)
}
\keyword{category}% at least one, from doc/KEYWORDS



