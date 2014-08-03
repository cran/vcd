\name{table2d_summary}
\alias{table2d_summary}
\alias{print.table2d_summary}
\title{Summary of a 2-way Table}
\description{
  Prints a 2-way contingency table along with percentages, marginal,
  and conditional distributions.
}
\usage{
table2d_summary(object, margins = TRUE, percentages = FALSE,
  conditionals = c("none", "row", "column"), \dots)
}
\arguments{
  \item{object}{a \eqn{r \times c}{r x c}-contingency table}
  \item{margins}{if \code{TRUE}, marginal distributions are computed.}
  \item{percentages}{if \code{TRUE}, relative frequencies are computed.}
  \item{conditionals}{if not \code{"none"}, the conditional
    distributions, given the row/column factor, are computed.}
  \item{\dots}{currently not used.}
}
\value{
  Returns invisibly a \eqn{r \times c \times k}{r x c x k}
  table, \eqn{k} depending on the amount of choices (at most 3).
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{mar_table}},
  \code{\link{prop.table}},
  \code{\link{independence_table}}
}
\examples{
data("UCBAdmissions")
table2d_summary(margin.table(UCBAdmissions, 1:2))
}
\keyword{category}
