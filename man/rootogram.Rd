\name{rootogram}
\alias{rootogram}
\alias{rootogram.default}
\alias{rootogram.goodfit}
\title{Rootograms}
\description{
Rootograms of observed and fitted values.
}
\usage{
\method{rootogram}{default}(x, fitted, names = NULL, scale = c("sqrt", "raw"),
                            type = c("hanging", "standing", "deviation"),
			    bar.col = grey(0.7), line.col = 2,
			    xlab = NULL, ylab = NULL, ylim = NULL, ...)
}
\arguments{
  \item{x}{either a vector or a 1-way table of frequencies.}
  \item{fitted}{a vector of fitted frequencies.}
  \item{names}{a vector of names passed to barplot, if set to
    \code{NULL} the names of \code{x} are used.}
  \item{scale}{a character string indicating wether the values should be
    plotted on the raw or square root scale.}
  \item{type}{a character string indicating if the bars for the observed
    frequencies should be hanging or standing or indicate the deviation
    between observed and fitted frequencies.}
  \item{bar.col}{bar color (observed frequencies).}
  \item{line.col}{line color (fitted frequencies).}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{ylim}{limits for the y axis.}
  \item{...}{further arguments passed to \code{\link{barplot}}.}
}
\details{
The observed frequencies are displayed as bars and the fitted frequencies as
a line. By default a log scale is used to make the smaller frequencies more
visible.
}
\author{
  Achim Zeileis
}
\references{
J. W. Tukey (1977), \emph{Exploratory Data Analysis}. Addison Wesley,
Reading, MA.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\examples{
## Simulated data examples:
dummy <- rnbinom(200, size = 1.5, prob = 0.8)
observed <- table(dummy)
fitted1 <- dnbinom(as.numeric(names(observed)), size = 1.5, prob = 0.8) * sum(observed)
fitted2 <- dnbinom(as.numeric(names(observed)), size = 2, prob = 0.6) * sum(observed)
rootogram(observed, fitted1)
rootogram(observed, fitted2)

## Real data examples:
data(HorseKicks)
HK.fit <- goodfit(HorseKicks)
summary(HK.fit)
plot(HK.fit)
## or equivalently
rootogram(HK.fit)

data(Federalist)
F.fit <- goodfit(Federalist, type = "nbinomial")
summary(F.fit)
plot(F.fit)
}
\keyword{hplot}
