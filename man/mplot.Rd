\name{mplot}
\alias{mplot}
\title{Multiple Grid plots}
\description{
  combines severals grid-based plots in a multi-panel-layout.
}
\usage{
mplot(..., .list = list(),
      layout = NULL, cex = NULL,
      main = NULL, gp_main = gpar(fontsize = 20),
      sub = NULL, gp_sub = gpar(fontsize = 15),
      keep_aspect_ratio = TRUE)
}
\arguments{
  \item{\dots, .list}{A list of objects inheriting from class \code{"grob"},
    or having a \code{"grob"} attribute containing such an object.}
  \item{layout}{integer vector of length 2 giving the number of rows and
    columns. If \code{NULL}, the values will be
    guessed using some heuristics from the number of objects supplied in
    \dots.}
  \item{cex}{Scaling factor for the fonts in the subplots. If
    \code{NULL}, the value is calculated as the inverse square root of
    the row number.}
  \item{main, sub}{Optional main and sub title, respectively.}
  \item{gp_main, gp_sub}{Optional objects of class \code{"gpar"} specifying
    the graphical parameters for the main and sub title, respectively.}
  \item{keep_aspect_ratio}{logical; should the aspect ratio of the plots
  be fixed?}
}
\value{
  None.
}
\details{
  This is a convenience function for producing multi-panel plots from
  grid-based displays, especially those produced by the vcd methods. The
  layout (number of rows and columns) is guessed from the amount of
  supplied objects, if not supplied. Currently, the vcd plotting
  functions do not return grob objects by default---this might change in
  the future. Also, some of them will return the grob object as a
  \code{"grob"} attribute, attached to the currently returned object.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
mplot(mosaic(Titanic, return_grob = TRUE),
      assoc(Titanic), return_grob = TRUE)

A = mosaic(Titanic, return_grob = TRUE)
B = mosaic(Titanic, type = "expected", return_grob = TRUE)
mplot(A, B)

mplot(sieve(SexualFun, return_grob = TRUE),
      agreementplot(SexualFun, return_grob = TRUE),
      main = "Sexual Fun")

mplot(A, grid.circle())
}