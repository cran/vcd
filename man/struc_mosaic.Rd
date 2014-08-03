\name{struc_mosaic}
\alias{struc_mosaic}
\title{Core-generating Function for Mosaic Plots}
\description{
  Core-generating function for \code{strucplot} returning a function
  producing mosaic plots.
}
\usage{
struc_mosaic(zero_size = 0.5, zero_split = FALSE, zero_shade = TRUE,
             zero_gp = gpar(col = 0), panel = NULL)
}
\arguments{
  \item{zero_size}{size of the bullets used for zero-entries in the
    contingency table (if 0, no bullets are drawn).} 
  \item{zero_split}{logical controlling whether zero cells should be
    further split. If \code{FALSE} and \code{zero_shade} is
    \code{FALSE}, only one bullet is drawn (centered) for unsplit
    zero cells. If \code{FALSE} and \code{zero_shade} is
    \code{TRUE}, a bullet for each zero cell is drawn to allow, e.g.,
    residual-based shadings to be effective also for zero cells.}
  \item{zero_shade}{logical controlling whether zero bullets should be
    shaded.}
  \item{zero_gp}{object of class \code{"gpar"} used
    for zero bullets in case they are \emph{not} shaded.}
  \item{panel}{Optional function with arguments: \code{residuals},
    \code{observed}, \code{expected}, \code{index}, \code{gp}, and
    \code{name} called by the
    \code{struc_mosaic} workhorse for each tile that is drawn in the
    mosaic. \code{index} is an integer vector with the tile's
    coordinates in the contingency table, \code{gp} a \code{gpar} object
    for the tile, and \code{name} a label to be assigned to the drawn
    grid object.}
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
  \url{http://datavis.ca/sugi/sugi17.pdf}

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
## mosaic plot with large zeros
strucplot(Titanic, core = struc_mosaic(zero_size = 1))
}
\keyword{hplot}
