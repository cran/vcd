\name{tile}
\alias{tile}
\alias{tile.default}
\alias{tile.formula}
\title{Tile Plot}
\description{
  Plots a tile display.
}
\usage{
\method{tile}{default}(x,
         tile_type = c("squaredarea", "area", "height", "width"),
         halign = c("left", "center", "right"),
         valign = c("bottom", "center", "top"),
         split_vertical = NULL,
         shade = FALSE,
         spacing = spacing_equal(unit(1, "lines")),
         set_labels = NULL,
         margins = unit(3, "lines"),
         keep_aspect_ratio = FALSE,
         legend_width = NULL,
         squared_tiles = TRUE,
         main = NULL, sub = NULL, ...)
\method{tile}{formula}(formula, data,
  \dots, main = NULL, sub = NULL, subset = NULL, na.action = NULL)
}
\arguments{
  \item{x}{a contingency table, or an object coercible to one.}
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}}.
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s. Ignored if \code{data} is a contingency
    table.}
  \item{tile_type}{character string indicating how the tiles should reflect
    the table frequencies (see details).}
  \item{halign, valign}{character string specifying the horizontal and
    vertical alignment of the tiles.}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where
    \eqn{k} is the number of margins of \code{x} (values are recycled as needed).
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits. Default is \code{FALSE.}}
  \item{spacing}{spacing object, spacing function, or corresponding
    generating function (see \code{\link{strucplot}} for more
    information).}
  \item{set_labels}{An optional character vector with named components
    replacing the so-specified variable names. The component names must
    exacty match the variable names to be replaced.}
  \item{shade}{logical specifying whether shading should be enabled or
    not (see \code{\link{strucplot}}).}
  \item{margins}{either an object of class \code{"unit"} of length 4, or
    a numeric vector of length 4. The elements are recycled as needed.
    The four components specify the top, right,
    bottom, and left margin of the plot, respectively.
    When a numeric vector is supplied, the numbers are interpreted as
    \code{"lines"} units.  In addition, the unit or numeric vector
    may have named arguments
    (\samp{top}, \samp{right}, \samp{bottom}, and \samp{left}), in which
    case the non-named arguments specify the default values (recycled as
    needed), overloaded by the named arguments.}
  \item{legend_width}{An object of class \code{"unit"} of length
    1 specifying the width of the legend (if any). Default: 5 lines.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be
    fixed or not. The default is \code{FALSE} to enable the creation of
    squared tiles.}
  \item{squared_tiles}{logical indicating whether white space should be
    added as needed to rows or columns to obtain squared tiles in case of an
    unequal number of row and column labels.}
  \item{main, sub}{either a logical, or a character string used for plotting
    the main (sub) title.  If logical and \code{TRUE}, the
    name of the \code{data} object is used.}
  \item{\dots}{Other arguments passed to \code{\link{strucplot}}}
}
\details{
  A tile plot is a matrix of tiles. For each tile, either the \code{"width"},
  \code{"height"}, \code{"area"}, or squared area is proportional to the
  corresponding entry. The first three options allow
  column-wise, row-wise and overall comparisons, respectively. The last
  variant allows to compare the tiles both column-wise and row-wise,
  considering either the width or the height, respectively.

  In contrast to other high-level strucplot functions, \code{tile}
  also accepts a table with duplicated levels (see examples). In this
  case, artificial dimnames will be created, and the actual ones are
  drawn using \code{set_labels}.
  
  Note that multiway-tables are first \dQuote{flattened} using
  \code{structable}.
}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\seealso{
  \code{\link{assoc}},
  \code{\link{strucplot}},
  \code{\link{mosaic}},
  \code{\link{structable}},
}
\examples{
data("Titanic")

## default plot
tile(Titanic)
tile(Titanic, type = "expected")
tile(Titanic, shade = TRUE)

## some variations
tile(Titanic, tile_type = "width", squared_tiles = FALSE)
tile(Titanic, tile_type = "height", squared_tiles = FALSE)
tile(Titanic, tile_type = "area", halign = "center", valign = "center")

## repeat levels
tile(Titanic[,,,c(1,2,1,2)])
}
\keyword{hplot}
