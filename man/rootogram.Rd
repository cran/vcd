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
  rect_gp = gpar(fill = "lightgray"), lines_gp = gpar(col = "red"),
  points_gp = gpar(col = "red"), pch = 19,
  xlab = NULL, ylab = NULL, ylim = NULL,
  name = "rootogram", newpage = TRUE, pop = TRUE, \dots)
}
\arguments{
  \item{x}{either a vector or a 1-way table of frequencies.}
  \item{fitted}{a vector of fitted frequencies.}
  \item{names}{a vector of names passed to \code{\link{grid_barplot}}, if set
    to \code{NULL} the names of \code{x} are used.}
  \item{scale}{a character string indicating whether the values should be
    plotted on the raw or square root scale.}
  \item{type}{a character string indicating if the bars for the observed
    frequencies should be hanging or standing or indicate the deviation
    between observed and fitted frequencies.}
  \item{rect_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the rectangles.}
  \item{lines_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the lines.}
  \item{points_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the points.}
  \item{pch}{plotting character for the points.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{ylim}{limits for the y axis.}
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
  \item{\dots}{further arguments passed to \code{\link{grid_barplot}}.}
}
\details{
  The observed frequencies are displayed as bars and the fitted
  frequencies as a line.  By default a sqrt scale is used to make the
  smaller frequencies more visible.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\references{
  J. W. Tukey (1977),
  \emph{Exploratory Data Analysis}.
  Addison Wesley, Reading, MA.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\seealso{\code{\link{grid_barplot}}}
\examples{
## Simulated data examples:
dummy <- rnbinom(200, size = 1.5, prob = 0.8)
observed <- table(dummy)
fitted1 <- dnbinom(as.numeric(names(observed)),
                   size = 1.5, prob = 0.8) * sum(observed)
fitted2 <- dnbinom(as.numeric(names(observed)),
                   size = 2, prob = 0.6) * sum(observed)
rootogram(observed, fitted1)
rootogram(observed, fitted2)

## Real data examples:
data("HorseKicks")
HK.fit <- goodfit(HorseKicks)
summary(HK.fit)
plot(HK.fit)
## or equivalently
rootogram(HK.fit)

data("Federalist")
F.fit <- goodfit(Federalist, type = "nbinomial")
summary(F.fit)
plot(F.fit)
}
\keyword{hplot}
