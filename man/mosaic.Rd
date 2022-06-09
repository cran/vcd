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
  spacing_args = list(), gp = NULL, expected = NULL, shade = NULL,
  highlighting = NULL, highlighting_fill = rev(gray.colors(tail(dim(x), 1))),
  highlighting_direction = NULL,
  zero_size = 0.5, zero_split = FALSE, zero_shade = NULL,
  zero_gp = gpar(col = 0), panel = NULL, main = NULL, sub = NULL, \dots)
\method{mosaic}{formula}(formula, data, highlighting = NULL,
  \dots, main = NULL, sub = NULL, subset = NULL, na.action = NULL)
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
    used first for splitting. If any, a specified response variable will be
    highlighted in the cells.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s. Ignored if \code{data} is a contingency table.}
  \item{zero_size}{size of the bullets used for zero entries (if 0, no
    bullets are drawn).}
  \item{zero_split}{logical controlling whether zero cells should be
    further split. If \code{FALSE} and \code{zero_shade} is
    \code{FALSE}, only one bullet is drawn (centered) for unsplit
    zero cells. If \code{FALSE} and \code{zero_shade} is
    \code{TRUE}, a bullet for each zero cell is drawn to allow, e.g.,
    residual-based shadings to be effective also for zero cells.}
  \item{zero_shade}{logical controlling whether zero bullets should be
    shaded. The default is \code{TRUE} if \code{shade} is \code{TRUE}
    or \code{expected} is not null or \code{gp} is not null, and
   \code{FALSE} otherwise.}
  \item{zero_gp}{object of class \code{"gpar"} used
    for zero bullets in case they are \emph{not} shaded.}
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
    \code{spacing_increase} for more dimensions, and
    \code{spacing_conditional} if conditioning variables are specified
    using \code{condvars} or the formula interface.}
  \item{spacing_args}{list of arguments for the generating function, if
    specified (see \code{\link{strucplot}} for more information).}
  \item{gp}{object of class \code{"gpar"}, shading function or a
    corresponding generating function (see details and
    \code{\link{shadings}}). Components of \code{"gpar"}
    objects are recycled as needed along the last splitting dimension.
    Ignored if \code{shade = FALSE}.}
  \item{shade}{logical specifying whether \code{gp} should be used or not
    (see \code{gp}). If \code{TRUE} and \code{expected} is unspecified,
    a default model is fitted: if \code{condvars} (see
    \code{\link{strucplot}}) is specified, a
    corresponding conditional independence model, and else the total
    independence model.}
  \item{expected}{optionally, an array of expected values of the same dimension
    as \code{x}, or alternatively the corresponding independence model specification
    as used by \code{\link[stats]{loglin}} or \code{\link[MASS]{loglm}}
    (see \code{\link{strucplot}}).}
  \item{highlighting}{character vector or integer specifying a variable
    to be highlighted in the cells.}
  \item{highlighting_fill}{color vector or palette function
    used for a highlighted variable, if any.}
  \item{highlighting_direction}{Either \code{"left"}, \code{"right"},
    \code{"top"}, or \code{"bottom"} specifying the direction of
    highlighting in the cells.}
  \item{panel}{Optional function with arguments: \code{residuals},
    \code{observed}, \code{expected}, \code{index}, \code{gp}, and
    \code{name} called by the
    \code{struc_mosaic} workhorse for each tile that is drawn in the
    mosaic. \code{index} is an integer vector with the tile's
    coordinates in the contingency table, \code{gp} a \code{gpar} object
    for the tile, and \code{name} a label to be assigned to the drawn
    grid object.}
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
  Hartigan, J.A., and Kleiner, B. (1984),
  A mosaic of television ratings.
  \emph{The American Statistician},
  \bold{38}, 32--35.

  Emerson, J. W. (1998),
  Mosaic displays in S-PLUS: A general implementation and a case study.
  \emph{Statistical Computing and Graphics Newsletter (ASA)},
  \bold{9}, 1, 17--23.

  Friendly, M. (1994),
  Mosaic displays for multi-way contingency tables.
  \emph{Journal of the American Statistical Association}, \bold{89},
  190--200.

  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  \doi{10.18637/jss.v017.i03} and available as 
  \code{vignette("strucplot", package = "vcd")}.

  The home page of Michael Friendly
  (\url{http://datavis.ca}) provides
  information on various aspects of graphical methods for analyzing
  categorical data, including mosaic plots. In particular, there are
  many materials for his course \dQuote{Visualizing Categorical Data with SAS and R}
  at \url{http://datavis.ca/courses/VCD/}.
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
library(MASS)
data("Titanic")
mosaic(Titanic)

## Formula interface for tabulated data plus shading and legend:
mosaic(~ Sex + Age + Survived, data = Titanic,
  main = "Survival on the Titanic", shade = TRUE, legend = TRUE)

data("HairEyeColor")
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
data("mtcars")
mosaic(~ gear + carb, data = mtcars, shade = TRUE)

data("PreSex")
mosaic(PreSex, condvars = c(1,4))
mosaic(~ ExtramaritalSex + PremaritalSex | MaritalStatus + Gender,
       data = PreSex)

## Highlighting:
mosaic(Survived ~ ., data = Titanic)

data("Arthritis")
mosaic(Improved ~ Treatment | Sex, data = Arthritis, zero_size = 0)
mosaic(Improved ~ Treatment | Sex, data = Arthritis, zero_size = 0,
       highlighting_direction = "right")
}
\keyword{hplot}
