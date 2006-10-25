\name{struc_sieve}
\alias{struc_sieve}
\title{Core-generating Function for Sieve Plots}
\description{
  Core-generating function for \code{strucplot} returning a function
  producing sieve plots.
}
\usage{
struc_sieve(sievetype = c("observed","expected"))
}
\arguments{
  \item{sievetype}{logical indicating whether rectangles should be filled
    according to \code{observed} or \code{expected} frequencies.}
}
\details{
  This function is usually called by \code{\link{strucplot}} (typically
  when called by \code{\link{sieve}}) and returns a function used by
  \code{\link{strucplot}} to produce sieve plots. 
}
\value{
  A function with arguments:
  \item{residuals}{table of residuals.}
  \item{observed}{table of observed values.}
  \item{expected}{not used by \code{struc_sieve}.}
  \item{spacing}{object of class \code{"unit"} specifying the space
    between the tiles.}
  \item{gp}{list of \code{gpar} objects used for the drawing the tiles.}
  \item{split_vertical}{vector of logicals indicating, for each
    dimension of the table, the split direction.}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{sieve}},
  \code{\link{strucplot}},
  \code{\link{structable}}
}
\references{
  Riedwyl, H., and Schüpbach, M. (1994),
  Parquet diagram to plot contingency tables.
  In F. Faulbaum (ed.),
  \emph{Softstat '93: Advances in Statistical Software}, 293--299.
  Gustav Fischer, New York.

  Friendly, M. (2000),
  Visualizing Categorical Data,
  SAS Institute, Cary, NC.

  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.

}
\examples{
## Titanic data
data("Titanic")
strucplot(Titanic, core = struc_sieve)
}
\keyword{hplot}
