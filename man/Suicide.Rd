\name{Suicide}
\alias{Suicide}
\docType{data}
\title{Suicide Rates in Germany}
\description{
  Data from Heuer (1979) on suicide rates in West Germany classified by
  age, sex, and method of suicide.
}
\usage{
data(Suicide)
}
\format{
  A data frame with 306 observations and 6 variables.
  \describe{
    \item{Freq}{frequency of suicides.}
    \item{sex}{factor indicating sex (male, female).}
    \item{method}{factor indicating method used.}
    \item{age}{age (rounded).}
    \item{age.group}{factor.  Age classified into 5 groups.}
    \item{method2}{factor indicating method used (same as \code{method}
      but some levels are merged).}
  }
}
\references{
  J. Heuer (1979),
  \emph{Selbstmord bei Kindern und Jugendlichen}.
  Ernst Klett Verlag, Stuttgart.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  Michael Friendly (2000),
  Visualizing Categorical Data:
  \url{http://hotspur.psych.yorku.ca/ftp/sas/vcd/catdata/suicide.sas}
}
\examples{
data(Suicide)
structable(~ sex + method2 + age.group, data = Suicide)
}
\keyword{datasets}
