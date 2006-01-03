\name{co_table}
\alias{co_table}
\title{Compute Conditional Tables}
\description{
  For a contingency table in array form, compute a list of conditional
  tables given some margins.
}
\usage{
co_table(x, margin, collapse = ".")
}
\arguments{
  \item{x}{a contingency table in array form.}
  \item{margin}{margin index(es) or corresponding name(s) of the
    conditioning variables.}
  \item{collapse}{character used when collapsing level names
    (if more than 1 \code{margin} is specified).}
}
\details{
  This is essentially an interface to \code{\link[base]{"["}}
  which is more convenient for arrays of arbitrary dimension.   
}
\value{
  A list of the resulting conditional tables.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
data("HairEyeColor")
co_table(HairEyeColor, 1)
co_table(HairEyeColor, c("Hair", "Eye"))
co_table(HairEyeColor, 1:2, collapse = "")
}
\keyword{array}
