\name{woolf.test}
\alias{woolf.test}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Woolf test}
\description{
Test for homogeneity on 2 x 2 x k - tables over strata (i.e., if the log
odds ratios are the same in all strata).
}
\usage{
woolf.test(x)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A 2x2xk-table.}
}
\value{
     A list with class `"htest"' containing the following components: 
\item{statistic}{the value the chi-squared test statistic.}
\item{parameter}{the degrees of freedom of the approximate chi-squared
          distribution of the test statistic.}
\item{p.value}{the p-value for the test.}
\item{method}{a character string indicating the type of test performed.}
\item{data.name}{a character string giving the name(s) of the data.}
\item{observed}{the observed counts.}
\item{expected}{the expected counts under the null hypothesis.}
}

\seealso{\code{\link{mantelhaen.test}}}

\examples{
data(CoalMiners)
woolf.test(CoalMiners)
}
\keyword{htest}% at least one, from doc/KEYWORDS

