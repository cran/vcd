\name{ternaryplot}
\alias{ternaryplot}
\title{Ternary Diagram}
\description{
  Visualizes compositional, 3-dimensional data in an equilateral
  triangle.
}
\usage{
ternaryplot(x, scale = 1, dimnames = NULL,
  dimnames_position = c("corner","edge","none"),
  dimnames_color = "black", id = NULL, id_color = "black",
  coordinates = FALSE, grid = TRUE, grid_color = "gray",
  labels = c("inside", "outside", "none"),
  labels_color = "darkgray", border = "black", bg = "white",
  pch = 19, cex = 1, prop_size = FALSE, col = "red",
  main = "ternary plot", newpage = TRUE, pop = TRUE, \dots)
}
\arguments{
  \item{x}{a matrix with three columns.}
  \item{scale}{row sums scale to be used.}
  \item{dimnames}{dimension labels (defaults to the column names of
    \code{x}).}
  \item{dimnames_position, dimnames_color}{position and color of
    dimension labels.}
  \item{id}{optional labels to be plotted below the plot
    symbols. \code{coordinates} and \code{id} are mutual exclusive.}
  \item{id_color}{color of these labels.}
  \item{coordinates}{if \code{TRUE}, the coordinates of the points are
    plotted below them. \code{coordinates} and \code{id} are mutual
    exclusive.}
  \item{grid}{if \code{TRUE}, a grid is plotted. May optionally
    be a string indicating the line type (default: \code{"dotted"}).}
  \item{grid_color}{grid color.}
  \item{labels, labels_color}{position and color of the grid labels.}
  \item{border}{color of the triangle border.}
  \item{bg}{triangle background.}
  \item{pch}{plotting character. Defaults to filled dots.}
  \item{cex}{a numerical value giving the amount by which plotting text
    and symbols should be scaled relative to the default.  Ignored for
    the symbol size if \code{prop_size} is not \code{FALSE}.}
  \item{prop_size}{if \code{TRUE}, the symbol size is plotted
    proportional to the row sum of the three variables, i.e., represents
    the weight of the observation.}
  \item{col}{plotting color.}
  \item{main}{main title.}
  \item{newpage}{if \code{TRUE}, the plot will appear on a new graphics
    page.} 
  \item{pop}{logical; if \code{TRUE}, all newly generated viewports are popped after plotting.}
  \item{\dots}{additional graphics parameters (see \code{par})}
}
\details{
  A points' coordinates are found by computing the gravity center of
  mass points using the data entries as weights.  Thus, the coordinates
  of a point \eqn{P(a,b,c)}, \eqn{a + b + c = 1}, are: \eqn{P(b + c/2, c
    \sqrt{3}/2)}{P(b + c/2, c * sqrt(3)/2)}.
}
\references{
  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\examples{
data(Arthritis)

## Build table by crossing Treatment and Sex
tab <- as.table(xtabs(~ I(Sex:Treatment) + Improved, data = Arthritis))

## Mark groups
col <- c("red", "red", "blue", "blue")
pch <- c(1, 19, 1, 19)

## plot
ternaryplot(
            tab,
            col = col,
            pch = pch,
            prop_size = TRUE,
            bg = "lightgray",
            grid_color = "white",
            labels_color = "white",
            main = "Arthritis Treatment Data"
           )
## legend
grid_legend(0.8, 0.7, pch, col, rownames(tab), title = "GROUP")

## Titanic
data(Lifeboats)
attach(Lifeboats)
ternaryplot(
  Lifeboats[,4:6],
  pch = ifelse(side == "Port", 1, 19),
  col = ifelse(side == "Port", "red", "blue"),
  id  = ifelse(men / total > 0.1, as.character(boat), NA),
  prop_size = 2,
  dimnames_position = "edge",
  main = "Lifeboats on Titanic"
)
grid_legend(0.8, 0.9, c(1, 19),
  c("red", "blue"), c("Port", "Starboard"),
  title = "SIDE")


## Hitters
data(Hitters)
attach(Hitters)
colors <- c("black","red","green","blue","red","black","blue")
pch <- substr(levels(Positions), 1, 1)
ternaryplot(
  Hitters[,2:4],
  pch = as.character(Positions),
  col = colors[as.numeric(Positions)],
  main = "Baseball Hitters Data"
)
grid_legend(0.8, 0.9, pch, colors, levels(Positions),
  title = "POSITION(S)")

}
\keyword{hplot}
