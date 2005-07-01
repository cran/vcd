\name{OvaryCancer}
\alias{OvaryCancer}
\docType{data}
\title{Ovary Cancer Data}
\description{
  Data from Obel (1975) about a retrospective study of ovary cancer
  carried out in 1973.  Information was obtained from 299 women, who
  were operated for ovary cancer 10 years before.
}
\usage{
data(OvaryCancer)
}
\format{
  A data frame with 16 observations and 5 variables.
  \describe{
    \item{Freq}{frequency.}
    \item{stage}{factor indicating the stage of the cancer at the
      time of operation (early, advanced).}
    \item{operation}{factor indicating type of operation (radical,
      limited).}
    \item{survival}{factor indicating survival status after 10 years
      (yes, no).}
    \item{xray}{factor indicating whether X-ray treatment was received
      (yes, no).}
  }
}

\references{
  E. B. Obel (1975),
  A Comparative Study of Patients with Cancer of the Ovary Who Have
  Survived More or Less Than 10 Years.
  \emph{Acta Obstetricia et Gynecologica Scandinavica}, 
  \bold{55}, 429-439.

  E. B. Andersen (1991),
  \emph{The Statistical Analysis of Categorical Data}.
  2nd edition.
  Springer-Verlag, Berlin.
}

\source{
  E. B. Andersen (1991),
  The Statistical Analysis of Categorical Data, Table 6.4.
}

\examples{
data(OvaryCancer)
structable(~ ., data = OvaryCancer)
}

\keyword{datasets}
