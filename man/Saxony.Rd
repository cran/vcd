\name{Saxony}
\alias{Saxony}
\docType{data}
\title{Families in Saxony}
\description{
  Data from Geissler, cited in Sokal & Rohlf (1969) and Lindsey (1995)
  on gender distributions in families in Saxony in the 19th century.
}
\usage{
data(Saxony)
}
\format{
  A 1-way table giving the number of male children in 6115 families of
  size 12.  The variable and its levels are
  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab nMales \tab 0, 1, \dots, 12 \cr
  }
}

\references{
  J. K. Lindsey (1995),
  \emph{Analysis of Frequency and Count Data}.
  Oxford University Press, Oxford, UK.

  R. R. Sokal & F. J. Rohlf (1969),
  \emph{Biometry. The Principles and Practice of Statistics}.
  W. H. Freeman, San Francisco, CA.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  M. Friendly (2000),
  Visualizing Categorical Data, pages 40--42.
}
\examples{
data(Saxony)
gf <- goodfit(Saxony, type = "binomial")
summary(gf)
plot(gf)
}
\keyword{datasets}
