\name{UKSoccer}
\alias{UKSoccer}
\docType{data}
\title{UK Soccer Scores}
\description{
Data from Lee (1997), on the goals scored by Home and Away teams
in the Premier Football Leage, 1995/6 season.
}

\usage{data(UKSoccer)}

\format{A 2-dimensional array resulting from cross-tabulating the number
  of goals scored in 380 games.  The variables and their levels are as
  follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab Home \tab 0, 1, \dots, 4 \cr
    2 \tab Away \tab 0, 1, \dots, 4
  }
}

\references{
A. J. Lee (1997), Modelling scores in the Premier League: Is Manchester
United really the best?, \emph{Chance}, \bold{10}(1), 15-19.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data, p. 27.
}

\seealso{\code{\link{Bundesliga}}}

\examples{
data(UKSoccer)

mosaicplot(UKSoccer, shade = TRUE, main = "UK Soccer Scores")
}

\keyword{datasets}
