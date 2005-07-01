\name{strucplot}
\alias{strucplot}
\title{Structured Displays of Contingency Tables}
\description{
  This modular function visualizes certain aspects of high-dimensional
  contingency tables in a hierarchical way.
}
\usage{
strucplot(x, residuals = NULL, expected = NULL,
  condvars = NULL, shade = NULL, type = c("observed", "expected"),
  residuals_type = c("Pearson", "deviance", "FT"), split_vertical = TRUE, 
  spacing = spacing_equal, spacing_args = list(),
  gp = NULL, gp_args = list(),
  labeling = labeling_text, labeling_args = list(),
  panel = struc_mosaic, panel_args = list(),
  legend = NULL, legend_args = list(),
  main = NULL, sub = NULL, margins = rep.int(2.5, 4),
  legend_width = unit(0.15, "npc"), title_gp = gpar(fontsize = 20),
  newpage = TRUE, pop = TRUE, keep_aspect_ratio = TRUE)
}
\arguments{
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames} attribute.}
  \item{residuals}{optionally, an array of residuals of the same dimension
    as \code{x} (see details).}
  \item{expected}{optionally, an array of expected values of the same dimension
    as \code{x}, or alternatively the corresponding independence model specification
    as used by \code{\link[stats]{loglin}} or \code{\link[MASS]{loglm}} (see details).}
  \item{condvars}{number of conditioning variables, if any; those are
    expected to be ordered first in the table.
    This information is used for computing the expected values, and is
    also passed to the spacing functions (see \code{\link{spacings}}).}
  \item{shade}{logical specifying whether gp should be used or not
    (see \code{gp}). If \code{TRUE} and \code{expected} is unspecified,
    a default model is fitted: if \code{condvars} is specified, a
    corresponding conditional independence model, and else the total
    independence model.}
  \item{residuals_type}{a character string indicating the type of
    residuals to be computed when none are supplied.
    Must be one of \code{"pearson"} (giving components of Pearson's
    chi-squared), \code{"deviance"} (giving components of the likelihood
    ratio chi-squared), or \code{"FT"} for the Freeman-Tukey residuals.
    The value of this argument can be abbreviated.}
  \item{type}{a character string indicating whether the observed or the
    expected values of the table should be visualized.}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where
    \eqn{k} is the number of margins of \code{x} (values are recycled as needed).
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits.}
  \item{spacing}{spacing object, spacing function, or a corresponding
    generating function (see details and \code{\link{spacings}}).}
  \item{spacing_args}{list of arguments for the spacing-generating
    function, if specified.}
  \item{gp}{object of class \code{"gpar"}, shading function or a
    corresponding generating function (see details and
    \code{\link{shadings}}).
    Ignored if \code{shade = FALSE}.} 
  \item{gp_args}{list of arguments for the shading-generating function,
    if specified.}
  \item{labeling}{either a labeling function, or a corresponding
    generating function (see details and \code{\link{labelings}}.}
  \item{labeling_args}{list of arguments for the labeling-generating
    function, if specified.}
  \item{panel}{either a panel function, or a corresponding generating
    function (see details).  Currently, generating functions for
    mosaic plots (\code{\link{struc_mosaic}}) and association plots
    (\code{\link{struc_assoc}}) are provided.}
  \item{panel_args}{list of arguments for the panel-generating function,
    if specified.}
  \item{legend}{either a legend-generating function, or a legend
    function (see details and \code{\link{legends}}), or a logical.
    If \code{legend} is \code{NULL} or \code{TRUE} and \code{gp} is a
    function, legend defaults to \code{\link{legend_resbased}}.}
  \item{legend_args}{list of arguments for the legend-generating
    function, if specified.}
  \item{main}{either a logical, or a character string used for plotting
    the main title.  If \code{main} is a logical and \code{TRUE}, the
    name of the object supplied as \code{x} is used.}
  \item{sub}{a character string used for plotting the subtitle.}
  \item{margins}{either an object of class \code{"unit"} of length 4, or
    a numeric vector of length 4. The elements are recycled as needed.
    The four components specify the \code{top}, \code{right},
    \code{bottom}, and \code{left} margin of the plot, respectively.
    When a numeric vector is supplied, the numbers are interpreted as
    \code{"lines"} units.  In addition, it may have named arguments
    (\samp{top}, \samp{right}, \samp{bottom}, and \samp{left}), in which
    case the non-named arguments specify the default values (recycled as
    needed), overloaded by the named arguments.}
  \item{legend_width}{either an object of class \code{"unit"} of length
    1, or a numeric value interpreted as proportion of the available
    horizontal space, specifying the width of the legend (if any).}
  \item{pop}{logical indicating whether the generated viewport tree
    should be removed at the end of the drawing or not.}
  \item{title_gp}{object of class \code{"gpar"} containing the graphical
    parameters used for the main and sub titles.}
  \item{newpage}{logical indicating whether a new page should be created
    for the plot or not.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be
    fixed or not.}
}
\details{
  This function---usually called by higher-level functions such as
  \code{\link{assoc}} and \code{\link{mosaic}}---generates conditioning
  plots of contingency tables. First, it sets up a set of viewports for
  main- and subtitles, legend, and the actual plot region.  Then,
  residuals are computed as needed from observed and expected
  frequencies, where the expected frequencies are optionally computed
  for a specified independence model.  Finally, the specified functions
  for spacing, gp, main plot, legend, and labeling are called to produce
  the plot.  The function invisibly returns the \code{"structable"} object
  visualized.

  Most elements of the plot, such as the panel function, the spacing
  between the tiles, the shading of the tiles, the labeling, and the
  legend, are modularized in functions and specified as parameters.  For
  each element \emph{foo} (= \code{spacing}, \code{labeling}, \code{panel},
  or \code{legend}), \code{strucplot} takes two arguments:
  \var{foo} and \var{foo\_args}, which can be used to specify the
  parameters in the following alternative ways:

  \enumerate{
    \item Passing a suitable function to \var{foo} which subsequently
    will be called from \code{strucplot} to compute shadings, labelings,
    etc.
    \item Passing a corresponding \emph{generating} function to \var{foo},
    along with parameters passed to \var{foo\_args}, that generates such a
    function. Generating functions must inherit from class \code{"panel_generator"}.
    \item Except for the shading functions (\var{shading\_bar}),
    passing \var{foo(foo\_args)} to the \var{foo} argument.
    \item For shadings and spacings, passing the final parameter object
    itself; see the corresponding help pages for more details on the
    data structures.
  }
  
  If legends are drawn, a \sQuote{cinemascope}-like layout is used for
  the plot to preserve the 1:1 aspect ratio.

  If \code{type = "expected"}, the expected values are passed to the
  \code{observed} argument of the panel function, and the observed
  values to the \code{expected} argument.

  Although the \code{gp} argument is typically used for shading, it can
  be used for arbitrary modifications of the tiles' graphics parameters
  (e.g., for highlighting particular cells, etc.).

}
\note{
  The created viewports, as well as the tiles and bullets, are named and
  thus can conveniently modified after a plot has been drawn (and
  \code{pop = FALSE}).
}
\value{
  Invisibly, an object of class \code{"structable"} corresponding to the
  plot.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{assoc}},
  \code{\link{mosaic}},
  \code{\link{struc_assoc}},
  \code{\link{struc_mosaic}},
  \code{\link{structable}},
  \code{\link{doubledecker}},
  \code{\link{labelings}},
  \code{\link{shadings}},
  \code{\link{legends}},
  \code{\link{spacings}}
}
\examples{
strucplot(Titanic)
strucplot(Titanic, panel = struc_assoc)
strucplot(Titanic, spacing = spacing_increase,
                   spacing_args = list(start = 0.5, rate = 1.5))
strucplot(Titanic, spacing = spacing_increase(start = 0.5, rate = 1.5))

## modify a tile's color
strucplot(Titanic, pop = FALSE)
grid.edit("rect..Class.1st..Sex.Male..Age.Adult..Survived.Yes",
          gp = gpar(fill = "red"))
}
\keyword{hplot}
