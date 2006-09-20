\name{sieve}
\alias{sieve}
\alias{sieve.default}
\alias{sieve.formula}
\title{Extended Sieve Plots}
\description{
  (Extended) sieve displays for n-way contingency tables: plots rectangles with
  areas proportional to the expected cell frequencies and filled with a
  number of squares equal to the observed frequencies.  Thus, the
  densities visualize the deviations of the observed from the expected
  values.
}
\usage{
\method{sieve}{default}(x, condvars = NULL, gp = NULL, shade = NULL,
  legend = FALSE, split_vertical = NULL, direction = NULL, spacing = NULL,
  spacing_args = list(), sievetype = c("observed","expected"),
  main = NULL, sub = NULL, \dots)
\method{sieve}{formula}(formula, data, \dots, main = NULL, sub = NULL, subset = NULL)
}
\arguments{
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute.}
  \item{condvars}{vector of integers or character strings indicating
    conditioning variables, if any. The table will be permuted to order
    them first.}
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}.  For convenience, conditioning
    formulas can be specified; the conditioning variables will then be
    used first for splitting.  Formulas for sieve displays (unlike
    those for doubledecker plots) have no response variable.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{shade}{logical specifying whether \code{gp} should be used or not
    (see \code{gp}). If \code{TRUE} and \code{expected} is unspecified,
    a default model is fitted: if \code{condvars} is specified, a
    corresponding conditional independence model, and else the total
    independence model. If \code{shade} is \code{NULL} (default),
    \code{gp} is used if specified.}
  \item{sievetype}{logical indicating whether rectangles should be filled
    according to \code{observed} or \code{expected} frequencies.}
  \item{gp}{object of class \code{"gpar"}, shading function or a
    corresponding generating function (see details of \code{\link{strucplot}} and
    \code{\link{shadings}}). Components of \code{"gpar"}
    objects are recycled as needed along the last splitting
    dimension. The default is a modified version of
    \code{\link{shading_Friendly}}:
    if \code{sievetype} is \code{"observed"}, cells with positive
    residuals are painted red, and cells with negative residuals
    blue. If \code{sievetype} is \code{"expected"}, the sieves' color is
    gray. Ignored if \code{shade = FALSE}.} 
  \item{legend}{either a legend-generating function, a legend
    function (see details of \code{\link{strucplot}} and
    \code{\link{legends}}), or a logical value.
    If \code{legend} is \code{NULL} or \code{TRUE} and \code{gp} is a
    function, legend defaults to \code{\link{legend_resbased}}.}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where \eqn{k}
    is the number of margins of \code{x} (default: \code{FALSE}).
    Values are recycled as needed. 
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits.
    Ignored if \code{direction} is not \code{NULL}.}
  \item{direction}{character vector of length \eqn{k}, where \eqn{k} is the
    number of margins of \code{x} (values are recycled as needed).
    For each component, a value of \code{"h"} indicates that the tile(s)
    of the corresponding dimension should be split horizontally, whereas
    \code{"v"} indicates vertical split(s).}
  \item{spacing}{spacing object, spacing function, or corresponding
    generating function (see \code{\link{strucplot}} for more
    information).
    The default is no spacing at all if \code{x} has two dimensions,
    and  \code{spacing_increase} for more dimensions.}
  \item{spacing_args}{list of arguments for the generating function, if
    specified (see \code{\link{strucplot}} for more information).}
  \item{main, sub}{either a logical, or a character string used for plotting
    the main (sub) title.  If logical and \code{TRUE}, the
    name of the \code{data} object is used.}
  \item{\dots}{Other arguments passed to \code{\link{strucplot}}}
}
\details{
  \code{sieve} is a generic function which currently has a default method and a
  formula interface.  Both are high-level interfaces to the
  \code{\link{strucplot}} function, and produce (extended) sieve
  displays.  Most of the functionality is described there, such as
  specification of the independence model, labeling, legend, spacing,
  shading, and other graphical parameters.

  The layout is very flexible: the specification of shading, labeling,
  spacing, and legend is modularized (see \code{\link{strucplot}} for
  details).
}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\note{To be faithful to the original definition by Riedwyl & Schüpbach,
   the default is to have no spacing between the tiles for two-way tables.}
\references{
  H. Riedwyl & M. Schüpbach (1994),
  Parquet diagram to plot contingency tables.
  In F. Faulbaum (ed.),
  \emph{Softstat '93: Advances in Statistical Software}, 293--299.
  Gustav Fischer, New York.

  M. Friendly (2000),
  Visualizing Categorical Data,
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{assoc}},
  \code{\link{strucplot}},
  \code{\link{mosaic}},
  \code{\link{structable}},
  \code{\link{doubledecker}}
}
\examples{
data("HairEyeColor")

## aggregate over 'sex':
(tab <- margin.table(HairEyeColor, c(2,1)))

## plot expected values:
sieve(tab, sievetype = "expected", shade = TRUE)

## plot observed table:
sieve(tab, shade = TRUE)

## plot complete diagram:
sieve(HairEyeColor, shade = TRUE)

## an example for the formula interface:
data("VisualAcuity")
sieve(Freq ~ right + left,  data = VisualAcuity)

## example with observed values in the cells:
sieve(Titanic, pop = FALSE, shade = TRUE)
labeling_cells(text = Titanic, gp_text = gpar(fontface = 2))(Titanic)
}
\keyword{hplot}
