\name{SpaceShuttle}
\alias{SpaceShuttle}
\docType{data}
\title{Space Shuttle O-ring Failures}
\description{
  Data from Dalal et al. (1989) about O-ring failures in the NASA space
  shuttle program.  The damage index comes from a discussion of the data
  by Tufte (1997).
}
\usage{
data(SpaceShuttle)
}
\format{
  A data frame with 24 observations and 6 variables.
  \describe{
    \item{FlightNumber}{Number of space shuttle flight.}
    \item{Temperature}{temperature during start (in degrees F).}
    \item{Pressure}{pressure.}
    \item{Fail}{did any O-ring failures occur? (no, yes).}
    \item{nFailures}{how many (of six) 0-rings failed?.}
    \item{Damage}{damage index.}
  }
}
\references{
  S. Dalal, E. B. Fowlkes, B. Hoadly (1989),
  Risk analysis of the space shuttle: Pre-Challenger prediction of
  failure,
  \emph{Journal of the American Statistical Association},
  \bold{84}, 945--957.

  E. R. Tufte (1997),
  \emph{Visual Explanations: Images and Quantities, Evidence and
    Narrative}.
  Graphics Press, Cheshire, CT. 

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  Michael Friendly (2000),
  Visualizing Categorical Data:
  \url{http://hotspur.psych.yorku.ca/ftp/sas/vcd/catdata/orings.sas}
}
\examples{
data(SpaceShuttle)
plot(nFailures/6 ~ Temperature, data = SpaceShuttle,
     xlim = c(30, 81), ylim = c(0,1),
     main = "NASA Space Shuttle O-Ring Failures",
     ylab = "Estimated failure probability",
     pch = 19, col = 4)
fm <- glm(cbind(nFailures, 6 - nFailures) ~ Temperature,
          data = SpaceShuttle,
          family = binomial)
lines(30 : 81,
      predict(fm, data.frame(Temperature = 30 : 81), type = "re"),
      lwd = 2)
abline(v = 31, lty = 3)
}

\keyword{datasets}
