\name{grid_barplot}
\alias{grid_barplot}
\title{Barplot}
\description{
  Barplots of 1-way tables in grid.
}
\usage{
grid_barplot(height, width = 0.8, offset = 0,
  names = NULL, xlim = NULL, ylim = NULL, xlab = "", ylab = "", main = "",
  gp = gpar(fill = "lightgray"), name = "grid_barplot",
  newpage = TRUE, pop = FALSE)

}
\arguments{
  \item{height}{either a vector or a 1-way table of frequencies.}
  \item{width}{width of the bars (recycled if needed to the number of bars).}
  \item{offset}{offset of the bars (recycled if needed to the number of bars).}
  \item{names}{a vector of names for the bars, if set
    to \code{NULL} the names of \code{height} are used.}
  \item{xlim}{limits for the x axis.}
  \item{ylim}{limits for the y axis.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{main}{a title for the plot.}
  \item{gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the rectangles.}
  \item{name}{name of the plotting viewport.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
}
\details{
  \code{grid_barplot} mimics (some of) the features of \code{\link{barplot}},
  but currently it only supports 1-way tables.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
grid_barplot(sample(1:6), names = letters[1:6])
}
\keyword{hplot}
