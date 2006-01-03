\name{Bundesliga}
\alias{Bundesliga}
\docType{data}
\title{Ergebnisse der Fussball-Bundesliga}
\description{
  Results from the first German soccer league in the years 1995/6 and
  2001/2.
}
\usage{
data("Bundesliga")
}
\format{A data frame with 612 observations and 6 variables.
  \describe{
    \item{HomeTeam}{character. Name of the home team.}
    \item{AwayTeam}{character. Name of the away team.}
    \item{HomeGoals}{number of goals scored by the home team.}
    \item{AwayGoals}{number of goals scored by the away team.}
    \item{Round}{day of the game.}
    \item{Year}{year in which the season started.}
  }
}
\details{
  The data for the season 1995/6 is taken from Knorr-Held (1999),
  the data for the season 2001/2 was collected by Achim Zeileis.
}
\references{
  Leonhard Knorr-Held (1999),
  Dynamic rating of sports teams.
  SFB 386 \dQuote{Statistical Analysis of Discrete Structures},
  Discussion paper \bold{98}.
}
\source{
  SFB 386 \dQuote{Statistical Analysis of Discrete Structures}
  \url{http://www.stat.uni-muenchen.de/service/datenarchiv/bundesliga/bundesliga.html}
}
\seealso{
  \code{\link{UKSoccer}}
}
\examples{
data("Bundesliga")

## independence of home and away goals
tore95 <- xtabs(~ HomeGoals + AwayGoals, data = Bundesliga,
                subset = Year == 1995)
tore95
mosaic(tore95, gp = shading_max)

## number of goals per game poisson distributed?
ntore1 <- xtabs(~ HomeGoals, data = Bundesliga, subset = Year == 1995)
ntore2 <- xtabs(~ AwayGoals, data = Bundesliga, subset = Year == 1995)
ntore3 <- table(apply(subset(Bundesliga, Year == 1995)[,3:4], 1, sum))

gf1 <- goodfit(ntore1)
gf2 <- goodfit(ntore2)
gf3 <- goodfit(ntore3)

summary(gf1)
summary(gf2)
summary(gf3)
plot(gf1)
plot(gf2)
plot(gf3)

Ord_plot(ntore1)
distplot(ntore1)
}

\keyword{datasets}
