\name{summary.table}
\alias{summary.table}
\alias{print.summary.table}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary of a m x n-table}
\description{
Prints a m x n-contingency table along with percentages, marginal and conditional distributions.
}
\usage{
summary.table(object, margins = TRUE, percentages = FALSE,
        conditionals = c("none", "row", "column"), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{a m x n-contingency table}
  \item{margins}{if \code{TRUE}, marginal distributions are computed.}
  \item{percentages}{if \code{TRUE}, relative frequencies are computed.}
  \item{conditionals}{if not \code{"none"}, the conditioinal
    distributions, given the row/column factor, are computed.}
  \item{\dots}{currently not used.}
}
\value{
Returns invisibly a m x n x o - table, o depending on the amount of
choices (max 3).
}
\author{
  David Meyer\cr
  \email{david.meyer@ci.tuwien.ac.at}
}
\seealso{\code{\link{mar.table}}, \code{\link{prop.table}}}

\examples{
data(UCBAdmissions)
summary(margin.table(UCBAdmissions, 1:2))
}
\keyword{category}% at least one, from doc/KEYWORDS
