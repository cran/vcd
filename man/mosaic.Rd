\name{mosaic}
\alias{mosaic}
\alias{mosaic.default}
\alias{mosaic.formula}
\title{Extended Mosaic Plots}
\description{
  Plots (extended) mosaic displays.
}
\usage{
\method{mosaic}{default}(x, condvars = NULL,
  split_vertical = NULL, direction = NULL, spacing = NULL,
  spacing_args = list(), zero_size = 0.5, main = NULL, sub = NULL, \dots)
\method{mosaic}{formula}(formula, data, \dots, main = NULL, sub = NULL, subset = NULL)
}
\arguments{
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute, or an object
    of class \code{"structable"}.}
  \item{condvars}{vector of integers or character strings indicating
    conditioning variables, if any. The table will be permuted to order
    them first.}
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}.  For convenience, conditioning
    formulas can be specified; the conditioning variables will then be
    used first for splitting.  Formulas for mosaic displays (unlike
    those for doubledecker plots) have no response variable.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{zero_size}{size of the bullets used for zero entries (if 0, no
    bullets are drawn).}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where \eqn{k}
    is the number of margins of \code{x} (default: \code{FALSE}).
    Values are recycled as needed. 
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits. Ignored if \code{direction} is not \code{NULL}.}
  \item{direction}{character vector of length \eqn{k}, where \eqn{k} is the
    number of margins of \code{x} (values are recycled as needed).
    For each component, a value of \code{"h"} indicates that the tile(s)
    of the corresponding dimension should be split horizontally, whereas
    \code{"v"} indicates vertical split(s).}
  \item{spacing}{spacing object, spacing function, or corresponding
    generating function (see \code{\link{strucplot}} for more
    information).
    The default is \code{spacing_equal} if \code{x} has two dimensions,
    and  \code{spacing_increase} for more dimensions.}
  \item{spacing_args}{list of arguments for the generating function, if
    specified (see \code{\link{strucplot}} for more information).}
  \item{main, sub}{either a logical, or a character string used for plotting
    the main (sub) title.  If logical and \code{TRUE}, the
    name of the \code{data} object is used.}
  \item{\dots}{Other arguments passed to \code{\link{strucplot}}}
}
\details{
  Mosaic displays have been suggested in the statistical literature
  by Hartigan and Kleiner (1984) and have been extended by Friendly
  (1994). \code{\link[graphics]{mosaicplot}} is a base graphics
  implementation and \code{mosaic} is a much more flexible and extensible
  grid implementation.

  \code{mosaic} is a generic function which currently has a default method and a
  formula interface.  Both are high-level interfaces to the
  \code{\link{strucplot}} function, and produce (extended) mosaic
  displays.  Most of the functionality is described there, such as
  specification of the independence model, labeling, legend, spacing,
  shading, and other graphical parameters.

  A mosaic plot is an area proportional visualization of a (possibly
  higher-dimensional) table of expected frequencies.  It is composed of
  tiles (corresponding to the cells) created by recursive vertical and
  horizontal splits of a square.  The area of each tile is proportional
  to the corresponding cell entry, \emph{given} the dimensions of
  previous splits.

  An \emph{extended} mosaic plot, in addition, visualizes the fit of a
  particular log-linear model. Typically, this is done by residual-based
  shadings where color and/or outline of the tiles visualize sign, size
  and possibly significance of the corresponding residual.
  
  The layout is very flexible: the specification of shading, labeling,
  spacing, and legend is modularized (see \code{\link{strucplot}} for
  details).

  In contrast to the \code{\link[graphics]{mosaicplot}} function in
  \pkg{graphics}, the splits start with the \emph{horizontal} direction
  by default to match the printed output of \code{\link{structable}}.
}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\references{
  Hartigan, J.A., and Kleiner, B. (1984)
  A mosaic of television ratings.
  \emph{The American Statistician},
  \bold{38}, 32--35.

  Emerson, J. W. (1998)
  Mosaic displays in S-PLUS: A general implementation and a case study.
  \emph{Statistical Computing and Graphics Newsletter (ASA)},
  \bold{9}, 1, 17--23.

  Friendly, M. (1994)
  Mosaic displays for multi-way contingency tables.
  \emph{Journal of the American Statistical Association}, \bold{89},
  190--200.

  The home page of Michael Friendly
  (\url{http://www.math.yorku.ca/SCS/friendly.html}) provides
  information on various aspects of graphical methods for analyzing
  categorical data, including mosaic plots.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{assoc}},
  \code{\link{strucplot}},
  \code{\link{mosaicplot}},
  \code{\link{structable}},
  \code{\link{doubledecker}}
}
\examples{
data(Titanic)
mosaic(Titanic)

## Formula interface for tabulated data plus shading and legend:
mosaic(~ Sex + Age + Survived, data = Titanic,
  main = "Survival on the Titanic", shade = TRUE, legend = TRUE)

data(HairEyeColor)
mosaic(HairEyeColor, shade = TRUE)
## Independence model of hair and eye color and sex.  Indicates that
## there are significantly more blue eyed blond females than expected
## in the case of independence (and too few brown eyed blond females).

mosaic(HairEyeColor, shade = TRUE, expected = list(c(1,2), 3))
## Model of joint independence of sex from hair and eye color.  Males
## are underrepresented among people with brown hair and eyes, and are
## overrepresented among people with brown hair and blue eyes, but not
## "significantly".

## Formula interface for raw data: visualize crosstabulation of numbers
## of gears and carburettors in Motor Trend car data.
data(mtcars)
mosaic(~ gear + carb, data = mtcars, shade = TRUE)

data(PreSex)
mosaic(PreSex, condvars = c(1,4))
mosaic(~ ExtramaritalSex + PremaritalSex | MaritalStatus + Gender,
       data = PreSex)
}
\keyword{hplot}
