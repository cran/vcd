\name{Arthritis}
\alias{Arthritis}
\docType{data}
\title{Arthritis Treatment Data}
\description{
  Data from Koch \& Edwards (1988) from a double-blind clinical trial
  investigating a new treatment for rheumatoid arthritis.
}
\usage{data("Arthritis")}
\format{
  A data frame with 84 observations and 5 variables.
  \describe{
    \item{ID}{patient ID.}
    \item{Treatment}{factor indicating treatment (Placebo, Treated).}
    \item{Sex}{factor indicating sex (Female, Male).}
    \item{Age}{age of patient.}
    \item{Improved}{ordered factor indicating treatment outcome (None,
      Some, Marked).}
  }
}

\references{
  G. Koch \& S. Edwards (1988),
  Clinical efficiency trials with categorical data.
  In K. E. Peace (ed.), \emph{Biopharmaceutical Statistics for Drug
    Development}, 403--451.
  Marcel Dekker, New York.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  Michael Friendly (2000),
  Visualizing Categorical Data:
  \url{http://euclid.psych.yorku.ca/ftp/sas/vcd/catdata/arthrit.sas}
}
\examples{
data("Arthritis")
art <- xtabs(~ Treatment + Improved, data = Arthritis, subset = Sex == "Female")
art

mosaic(art, gp = shading_Friendly)
mosaic(art, gp = shading_max)
}

\keyword{datasets}
