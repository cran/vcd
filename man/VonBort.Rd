\name{VonBort}
\alias{VonBort}
\docType{data}
\title{Von Bortkiewicz Horse Kicks Data}
\description{
  Data from von Bortkiewicz (1898), given by Andrews \& Herzberg (1985),
  on number of deaths by horse or mule kicks in 14 corps of the Prussian
  army.
}
\usage{
data("VonBort")
}
\format{
  A data frame with 280 observations and 4 variables.
  \describe{
    \item{deaths}{number of deaths.}
    \item{year}{year of the deaths.}
    \item{corps}{factor indicating the corps.}
    \item{fisher}{factor indicating whether the corresponding corps was
      considered by Fisher (1925) or not.}
  }
}
\references{
  D. F. Andrews \& A. M. Herzberg (1985),
  \emph{Data: A Collection of Problems from Many Fields for the Student
    and Research Worker}.
  Springer-Verlag, New York, NY.

  R. A. Fisher (1925),
  \emph{Statistical Methods for Research Workers}.
  Oliver & Boyd, London.

  L. von Bortkiewicz (1898),
  \emph{Das Gesetz der kleinen Zahlen}.
  Teubner, Leipzig.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  Michael Friendly (2000),
  Visualizing Categorical Data:
  \url{http://euclid.psych.yorku.ca/ftp/sas/vcd/catdata/vonbort.sas}
}
\seealso{
  \code{\link{HorseKicks}} for a popular subsample.
}
\examples{
data("VonBort")
## HorseKicks data
xtabs(~ deaths, data = VonBort, subset = fisher == "yes")
}

\keyword{datasets}
