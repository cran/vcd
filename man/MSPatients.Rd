\name{MSPatients}
\alias{MSPatients}
\docType{data}
\title{Diagnosis of Multiple Sclerosis}
\description{
Data from Westlund & Kurland (1953) on the diagnosis of multiple sclerosis (MS):
two samples of patients, one from Winnipeg and one from New Orleans, were
each rated by two neurologists (one from each city) in four diagnostic
categories.
}

\usage{data(MSPatients)}

\format{A 3-dimensional array resulting from cross-tabulating 218
  observations on 3 variables.  The variables and their levels are as
  follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab New Orleans Neurologist \tab Certain, Probable, Possible, Doubtful \cr
    2 \tab Winnipeg Neurologist \tab Certain, Probable, Possible, Doubtful \cr
    3 \tab Patients \tab Winnipeg, New Orleans
  }
}

\references{
K. B. Westlund & L. T. Kurland (1953), Studies on multiple sclerosis in
Winnipeg, Manitoba and New Orleans, Lousiana, \emph{American Journal of
Hygiene}, \bold{57}, 380-396.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data:
\url{http://hotspur.psych.yorku.ca/ftp/sas/vcd/catdata/msdiag.sas}}

\examples{
data(MSPatients)

par(mfrow = c(1,2))
agreementplot(t(MSPatients[,,1]), main = "Winnipeg Patients")
agreementplot(t(MSPatients[,,2]), main = "New Orleans Patients")
par(mfrow = c(1,1))
}

\keyword{datasets}
