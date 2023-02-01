\name{JobSatisfaction}
\alias{JobSatisfaction}
\docType{data}
\title{Job Satisfaction Data}
\description{
  Data from Petersen (1968) about the job satisfaction of 715 blue
  collar workers, selected from Danish Industry in 1968.
}
\usage{
data("JobSatisfaction")
}
\format{
  A data frame with 8 observations and 4 variables.
  \describe{
    \item{Freq}{frequency.}
    \item{management}{factor indicating quality of management (bad,
      good).} 
    \item{supervisor}{factor indicating supervisor's job satisfaction
      (low, high).}
    \item{own}{factor indicating worker's own job satisfaction (low,
      high).}
  }
}

\references{
  E. B. Andersen (1991),
  \emph{The Statistical Analysis of Categorical Data}.
  2nd edition.
  Springer-Verlag, Berlin.

  E. Petersen (1968),
  \emph{Job Satisfaction in Denmark}.
  (In Danish).
  Mentalhygiejnisk Forlag, Copenhagen.
}

\source{
  E. B. Andersen (1991),
  The Statistical Analysis of Categorical Data, Table 5.4.
}

\examples{
data("JobSatisfaction")
structable(~ ., data = JobSatisfaction)

mantelhaen.test(xtabs(Freq ~ own + supervisor + management,
                      data = JobSatisfaction))

}

\keyword{datasets}
