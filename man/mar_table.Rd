\name{mar_table}
\alias{mar_table}
\title{Table with Marginal Sums}
\description{
  Adds row and column sums to a two-way table.
}
\usage{
mar_table(x)
}
\arguments{
  \item{x}{a two-way table.}
}
\value{
  A table with row and column totals added.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data("SexualFun")
mar_table(SexualFun)
}
\keyword{category}
