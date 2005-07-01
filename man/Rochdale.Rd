\name{Rochdale}
\alias{Rochdale}
\docType{data}
\title{Rochdale Data}
\description{
  Information on 665 households of Rochdale, Lancashire, UK.
  The study was conducted to identify influence factors on economical
  activity of wives. 
}
\usage{
data(Rochdale)
}
\format{
  A 8-dimensional array resulting from cross-tabulating 665 observations
  on 8 variables.  The variables and their levels are as follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab EconActive \tab yes, no \cr
    2 \tab Age \tab <38, >38 \cr
    3 \tab HusbandEmployed \tab yes, no \cr
    4 \tab Child \tab yes, no \cr
    5 \tab Education \tab yes, no \cr
    6 \tab HusbandEducation \tab yes, no \cr
    7 \tab Asian \tab yes, no \cr
    8 \tab HouseholdWorking \tab yes, no \cr
  }
}

\note{
  Many observations are missing: only 91 out of all 256 combinations
  contain information.
}
\source{
  Whittaker (1990).
}
\references{
  H. Hofmann (2003).
  Constructing and reading mosaicplots.
  \emph{Computational Statistics \& Data Analysis},
  \bold(43), 4, 565--580.

  J. Whittaker (1990),
  \emph{Graphical Models on Applied Multivariate Statistics},
  Wiley, New York.
}
\examples{
data(Rochdale)
mosaic(Rochdale)
}
\keyword{datasets}
