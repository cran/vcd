\name{Federalist}
\alias{Federalist}
\docType{data}
\title{`May' in Federalist Papers}
\description{
  Data from Mosteller & Wallace (1984) investigating the use of certain
  keywords (\sQuote{may} in this data set) to identify the author of 12
  disputed \sQuote{Federalist Papers} by Alexander Hamilton, John Jay
  and James Madison.
}
\usage{
data(Federalist)
}
\format{
  A 1-way table giving the number of occurrences of \sQuote{may} in 262
  blocks of text.  The variable and its levels are
  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab nMay \tab 0, 1, \dots, 6 \cr
  }
}
\references{
  F. Mosteller & D. L. Wallace (1984),
  \emph{Applied Bayesian and Classical Inference: The Case of the
    Federalist Papers}.
  Springer-Verlag, New York, NY. 

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\source{
  Michael Friendly (2000),
  Visualizing Categorical Data, page 19.
}
\examples{
data(Federalist)
gf <- goodfit(Federalist, type = "nbinomial")
summary(gf)
plot(gf)
}
\keyword{datasets}
