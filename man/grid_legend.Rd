\name{grid_legend}
\alias{grid_legend}
\title{Legend Function for grid Graphics}
\description{
  This function can be used to add legends to \emph{grid-based} plots.
}
\usage{
grid_legend(x, y, pch, col, labels, frame = TRUE, hgap = unit(0.5, "lines"),
  vgap = unit(0.3, "lines"), default_units = "lines", gp = gpar(),
  draw = TRUE, title = "Legend:")
}
\arguments{
  \item{x, y}{coordinates of the legend}
  \item{pch}{integer vector of plotting symbols}
  \item{col}{character vector of colors for the symbols}
  \item{labels}{character vector of labels corresponding to the symbols}
  \item{frame}{logical indicating whether the legend should have a
    border or not.}
  \item{hgap}{object of class \code{"unit"} specifying the space between
  symbols and labels}
  \item{vgap}{object of class \code{"unit"} specifying the space between
  the lines}
  \item{default_units}{character string indicating the default unit}
  \item{gp}{object of class \code{"gpar"} used for the legend}
  \item{draw}{logical indicating whether the legend be drawn or not.}
  \item{title}{character string indicating the plot's title}
}
\value{
  Invisibly, the legend as a \code{"grob"} object.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link[graphics]{legend}}
}
\examples{
data(Lifeboats)
attach(Lifeboats)
ternaryplot(Lifeboats[,4:6],
  pch = ifelse(side == "Port", 1, 19),
  col = ifelse(side == "Port", "red", "blue"),
  id  = ifelse(men / total > 0.1, as.character(boat), NA),
  prop_size = 2,
  dimnames_position = "edge",
  main = "Lifeboats on Titanic")
grid_legend(0.8, 0.9, c(1, 19), c("red", "blue"),
  c("Port", "Starboard"), title = "SIDE")
}
\keyword{hplot}
