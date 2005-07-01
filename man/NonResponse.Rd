\name{NonResponse}
\alias{NonResponse}
\docType{data}
\title{Non-Response Survey Data}
\description{
  Data about non-response for a Danish survey in 1965.
}
\usage{
data(NonResponse)
}
\format{
  A data frame with 12 observations and 4 variables.
  \describe{
    \item{Freq}{frequency.}
    \item{residence}{factor indicating whether residence was in
      Copenhagen, in a city outside Copenhagen or at the countryside
      (Copenhagen, City, Country).}
    \item{response}{factor indicating whether a response was given (yes,
      no).}
    \item{gender}{factor indicating gender (male, female).}
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
  The Statistical Analysis of Categorical Data, Table 5.17.
}
\examples{
data(NonResponse)
structable(~ ., data = NonResponse)
}
\keyword{datasets}
