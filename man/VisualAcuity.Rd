\name{VisualAcuity}
\alias{VisualAcuity}
\docType{data}
\title{Visual Acuity in Left and Right Eyes}
\description{
  Data from Kendall & Stuart (1961) on unaided vision among 3,242 men
  and 7,477 women, all aged 30-39 and employed in the U.K. Royal
  Ordnance factories 1943-1946.
}
\usage{
data("VisualAcuity")
}
\format{
  A data frame with 32 observations and 4 variables.
  \describe{
    \item{Freq}{frequency of visual acuity measurements.}
    \item{right}{visual acuity on right eye.}
    \item{left}{visual acuity on left eye.}
    \item{gender}{factor indicating gender of patient.}
  }
}
\references{
  M. G. Kendall & A. Stuart (1961),
  \emph{The Advanced Theory of Statistics}, Vol. 2.
  Griffin, London.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  M. Friendly (2000),
  Visualizing Categorical Data:
  \url{http://euclid.psych.yorku.ca/ftp/sas/vcd/catdata/vision.sas}
}
\examples{
data("VisualAcuity")
structable(~ gender + left + right, data = VisualAcuity)
sieve(Freq ~ left + right | gender, data = VisualAcuity, shade = TRUE)
cotabplot(Freq ~ left + right | gender, data = VisualAcuity,
          panel = cotab_agreementplot)
}

\keyword{datasets}
