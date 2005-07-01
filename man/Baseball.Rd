\name{Baseball}
\alias{Baseball}
\docType{data}
\title{Baseball Data}
\description{
  Baseball data.
}
\usage{
data(Baseball)
}
\format{
  A data frame with 322 observations and 25 variables.
  \describe{
    \item{name1}{player's first name.}
    \item{name2}{player's last name.}
    \item{atbat86}{times at Bat: number of official plate appearances by
      a hitter.  It counts as an official at-bat as long as the batter
      does not walk, sacrifice, get hit by a pitch or reach base due to
      catcher's interference.}
    \item{hits86}{hits.}
    \item{homer86}{home runs.}
    \item{runs86}{the number of runs scored by a player.  A run is
      scored by an offensive player who advances from batter to runner
      and touches first, second, third and home base in that order
      without being put out.} 
    \item{rbi86}{Runs Batted In: A hitter earns a run batted in when he
      drives in a run via a hit, walk, sacrifice (bunt or fly) fielder's
      choice, hit-batsman or on an error (when the official scorer rules
      that the run would have scored anyway).}
    \item{walks86}{A \dQuote{walk} (or \dQuote{base on balls}) is an
      award of first base granted to a batter who receives four pitches
      outside the strike zone.}
    \item{years}{Years in the Major Leagues.  Seems to count all years a
      player has actually played in the Major Leagues, not necessarily
      consecutive.}
    \item{atbat}{career times at bat.}
    \item{hits}{career hits.}
    \item{homeruns}{career homeruns.}
    \item{runs}{career runs.}
    \item{rbi}{career runs batted in.}
    \item{walks}{career walks.}
    \item{league86}{player's league.}
    \item{div86}{player's division.}
    \item{team86}{player's team.}
    \item{posit86}{player's position (see \code{\link{Hitters}}).}
    \item{outs86}{number of putouts (see \code{\link{Hitters}})}
    \item{assist86}{number of assists (see \code{\link{Hitters}})}
    \item{error86}{number of assists (see \code{\link{Hitters}})}
    \item{sal87}{annary salary on opening day (in USD 1000).}
    \item{league87}{league in 1987.}
    \item{team87}{team in 1987.}
  }
}

\references{
  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  SAS System for Statistical Graphics, First Edition, page A2.3
}

\seealso{\code{\link{Hitters}}}

\examples{
data(Baseball)
}
\keyword{datasets}
