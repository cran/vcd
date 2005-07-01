\name{RepVict}
\alias{RepVict}
\docType{data}
\title{Repeat Victimization Data}
\description{
  Data from Reiss (1980) given by Fienberg (1980) about instances of
  repeat victimization for households in the U.S. National Crime
  Survey.
}
\usage{
data(RepVict)
}
\format{
  A 2-dimensional array resulting from cross-tabulating victimization.
  The variables and their levels are as follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab First Victimization \tab Rape, Assault, Robbery, Pickpocket,
    Personal Larceny, \cr
    \tab \tab Burglary, Household Larceny, Auto Theft \cr
    2 \tab Second Victimization \tab Rape, Assault, Robbery, Pickpocket,
    Personal Larceny,\cr
    \tab \tab Burglary, Household Larceny, Auto Theft
  }
}

\references{
  S. E. Fienberg (1980),
  \emph{The Analysis of Cross-Classified Categorical Data},
  MIT Press, Cambridge, 2nd edition.

  A. J. J. Reiss (1980),
  Victim proneness by type of crime in repeat victimization.
  In S. E. Fienberg & A. J. J. Reiss (eds.),
  \emph{Indicators of Crime and Criminal Justice}.
  U.S. Government Printing Office, Washington, DC.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  Michael Friendly (2000),
  Visualizing Categorical Data, page 113.
}
\examples{
data(RepVict)

mosaic(RepVict[-c(4,7),-c(4,7)], gp = shading_max,
       main = "Repeat Victimization Data")
}
\keyword{datasets}
