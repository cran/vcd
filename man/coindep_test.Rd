\name{coindep_test}
\alias{coindep_test}
\alias{fitted.coindep_test}
\title{Test for (Conditional) Independence}
\description{
  Performs a test of (conditional) independence of 2 margins
  in a contingency table by simulation from the marginal
  distribution of the input table under (conditional)
  independence.
}
\usage{
coindep_test(x, margin = NULL, n = 1000, 
  indepfun = function(x) max(abs(x)), aggfun = max,
  alternative = c("greater", "less"),
  pearson = TRUE)
}
\arguments{
  \item{x}{a contingency table.}
  \item{margin}{margin index(es) or corresponding name(s) of the
    conditioning variables. Each resulting conditional table has to
    be a 2-way table.}
  \item{n}{number of (conditional) independence tables to be drawn.}
  \item{indepfun}{aggregation function capturing independence in
    (each conditional) 2-way table.}
  \item{aggfun}{aggregation function aggregating the test statistics
    computed by \code{indepfun}.}
  \item{alternative}{a character string specifying the alternative
    hypothesis; must be either \code{"greater"} (default) or
    \code{"less"} (and may be abbreviated.)}
  \item{pearson}{logical. Should the table of Pearson residuals under
    independence be computed and passed to \code{indepfun} (default)
    or the raw table of observed frequencies?}
}
\details{
  If \code{margin} is \code{NULL} this computes a simple independence
  statistic in a 2-way table. Alternatively, \code{margin} can give 
  several conditioning variables and then conditional independence in
  the resulting conditional table is tested.
  
  By default, this uses a (double) maximum statistic of Pearson residuals.
  By changing \code{indepfun} or \code{aggfun} a (maximum of) Pearson Chi-squared
  statistic(s) can be computed or just the usual Pearson Chi-squared statistics
  and so on. Other statistics can be computed by changing \code{pearson}
  to \code{FALSE}.
  
  The function uses \code{\link{r2dtable}} to simulate the distribution
  of the test statistic under the null.
}
\value{
  A list of class \code{"coindep_test"} inheriting from \code{"htest"}
  with following components: 
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the \eqn{p} value for the test.}
  \item{method}{a character string indicating the type of the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{observed}{observed table of frequencies}
  \item{expctd}{expected table of frequencies}
  \item{residuals}{corresponding Pearson residuals}
  \item{margin}{the \code{margin} used}
  \item{dist}{a vector of size \code{n} with simulated values of the
    distribution of the statistic under the null.}
  \item{qdist}{the corresponding quantile function (for computing 
    critical values).}
  \item{pdist}{the corresponding distribution function (for computing
    \eqn{p} values).}
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link{chisq.test}},
  \code{\link{fisher.test}},
  \code{\link{r2dtable}}
}

\examples{
TeaTasting <- matrix(c(3, 1, 1, 3), nr = 2,
                     dimnames = list(Guess = c("Milk", "Tea"),
                                     Truth = c("Milk", "Tea"))
)
## compute maximum statistic
coindep_test(TeaTasting)
## compute Chi-squared statistic
coindep_test(TeaTasting, indepfun = function(x) sum(x^2))
## use unconditional asymptotic distribution
chisq.test(TeaTasting, correct = FALSE)
chisq.test(TeaTasting)


data("UCBAdmissions")
## double maximum statistic
coindep_test(UCBAdmissions, margin = "Dept")
## maximum of Chi-squared statistics
coindep_test(UCBAdmissions, margin = "Dept", indepfun = function(x) sum(x^2))
## Pearson Chi-squared statistic
coindep_test(UCBAdmissions, margin = "Dept", indepfun = function(x) sum(x^2), aggfun = sum)
## use unconditional asymptotic distribution
loglm(~ Dept * (Gender + Admit), data = UCBAdmissions)
}
\keyword{htest}
