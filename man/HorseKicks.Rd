\name{HorseKicks}
\alias{HorseKicks}
\docType{data}
\title{Death by Horse Kicks}
\description{
  Data from von Bortkiewicz (1898), given by Andrews \& Herzberg (1985),
  on number of deaths by horse or mule kicks in 10 (of 14 reported)
  corps of the Prussian army.  4 corps were not considered by Fisher
  (1925) as they had a different organization. This data set is a
  popular subset of the \code{\link{VonBort}} data.
}
\usage{
data("HorseKicks")
}
\format{
  A 1-way table giving the number of deaths in 200 corps-years.  The
  variable and its levels are
  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab nDeaths \tab 0, 1, \dots, 4 \cr
  }
}

\references{
  D. F. Andrews & A. M. Herzberg (1985),
  \emph{Data: A Collection of Problems from Many Fields for the Student
    and Research Worker}.
  Springer-Verlag, New York, NY.

  R. A. Fisher (1925),
  \emph{Statistical Methods for Research Workers}.
  Oliver \& Boyd, London.

  L. von Bortkiewicz (1898),
  \emph{Das Gesetz der kleinen Zahlen}.
  Teubner, Leipzig.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  Michael Friendly (2000),
  Visualizing Categorical Data, page 18.
}
\seealso{
  \code{\link{VonBort}}
}
\examples{
data("HorseKicks")
gf <- goodfit(HorseKicks)
summary(gf)
plot(gf)
}

\keyword{datasets}
