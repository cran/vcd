\name{BrokenMarriage}
\alias{BrokenMarriage}
\docType{data}
\title{Broken Marriage Data}
\description{
  Data from the Danish Welfare Study about broken marriages or permanent
  relationships depending on gender and social rank.
}
\usage{
data("BrokenMarriage")
}
\format{
  A data frame with 20 observations and 4 variables.
  \describe{
    \item{Freq}{frequency.}
    \item{gender}{factor indicating gender (male, female).}
    \item{rank}{factor indicating social rank (I, II, III, IV, V).}
    \item{broken}{factor indicating whether the marriage or permanent
      relationship was broken (yes, no).}
  }
}
\references{
  E. B. Andersen (1991),
  \emph{The Statistical Analysis of Categorical Data}.
  2nd edition.
  Springer-Verlag, Berlin.
}
\source{
  E. B. Andersen (1991),
  The Statistical Analysis of Categorical Data, page 177.
}
\examples{
data("BrokenMarriage")
structable(~ ., data = BrokenMarriage)
}

\keyword{datasets}
