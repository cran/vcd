\name{struc_mosaic}
\alias{struc_mosaic}
\title{Panel-generating Function for Mosaic Plots}
\description{
  Panel-generating function for \code{strucplot} returning a function
  producing mosaic plots.
}
\usage{
struc_mosaic(zero_size = 0.5)
}
\arguments{
  \item{zero_size}{size of the bullets used for zero-entries in the
    contingency table (if 0, no bullets are drawn).} 
}
\details{
  This function is usually called by \code{\link{strucplot}} (typically
  when called by \code{\link{mosaic}}) and returns a function used by
  \code{\link{strucplot}} to produce mosaic plots. 
}
\value{
  A function with arguments:
  \item{residuals}{table of residuals.}
  \item{observed}{table of observed values.}
  \item{expected}{not used by \code{struc_mosaic}.}
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
  \code{\link{mosaic}},
  \code{\link{strucplot}},
  \code{\link{structable}}
}
\references{
  Cohen, A. (1980),
  On the graphical display of the significant components in a two-way
  contingency table.
  \emph{Communications in Statistics---Theory and Methods}, \bold{A9},
  1025--1041.

  Friendly, M. (1992),
  Graphical methods for categorical data.
  \emph{SAS User Group International Conference Proceedings}, \bold{17},
  190--200.
  \url{http://www.math.yorku.ca/SCS/sugi/sugi17-paper.html}
}
\example{
## Titanic data
data(Titanic)
## mosaic plot with large zeros
strucplot(Titanic, panel = struc_mosaic(zero_size = 1))
}
\keyword{hplot}
