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
  residuals_type = NULL, split_vertical = NULL, 
  spacing = spacing_equal, spacing_args = list(),
  gp = NULL, gp_args = list(),
  labeling = labeling_border, labeling_args = list(),
  core = struc_mosaic, core_args = list(),
  legend = NULL, legend_args = list(),
  main = NULL, sub = NULL, margins = unit(3, "lines"),
  title_margins = NULL, legend_width = NULL,
  main_gp = gpar(fontsize = 20), sub_gp = gpar(fontsize = 15),
  newpage = TRUE, pop = TRUE, keep_aspect_ratio = NULL, prefix = "", \dots)
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
  \item{shade}{logical specifying whether \code{gp} should be used or not
    (see \code{gp}). If \code{TRUE} and \code{expected} is unspecified,
    a default model is fitted: if \code{condvars} is specified, a
    corresponding conditional independence model, and else the total
    independence model.}
  \item{residuals_type}{a character string indicating the type of
    residuals to be computed when none are supplied.
    If \code{residuals} is \code{NULL}, \code{residuals_type} must
    be one of \code{"pearson"} (default; giving components of Pearson's
    chi-squared), \code{"deviance"} (giving components of the likelihood
    ratio chi-squared), or \code{"FT"} for the Freeman-Tukey residuals.
    The value of this argument can be abbreviated. If \code{residuals}
    are specified, the value of \code{residuals_type} is just passed
    \dQuote{as is} to the legend function.}
  \item{type}{a character string indicating whether the observed or the
    expected values of the table should be visualized.}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where
    \eqn{k} is the number of margins of \code{x} (values are recycled as needed).
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits. Default is \code{FALSE.}}
  \item{spacing}{spacing object, spacing function, or a corresponding
    generating function (see details and \code{\link{spacings}}).}
  \item{spacing_args}{list of arguments for the spacing-generating
    function, if specified.}
  \item{gp}{object of class \code{"gpar"}, shading function or a
    corresponding generating function (see details and
    \code{\link{shadings}}). Components of \code{"gpar"}
    objects are recycled as needed along the last splitting dimension.
    Ignored if \code{shade = FALSE}.} 
  \item{gp_args}{list of arguments for the shading-generating function,
    if specified.}
  \item{labeling}{either a logical, or a labeling function, or a corresponding
    generating function (see details and \code{\link{labelings}}. If
    \code{FALSE} or \code{NULL}, no labeling is produced.}
  \item{labeling_args}{list of arguments for the labeling-generating
    function, if specified.}
  \item{core}{either a core function, or a corresponding generating
    function (see details).  Currently, generating functions for
    mosaic plots (\code{\link{struc_mosaic}}), association plots
    (\code{\link{struc_assoc}}), and sieve plots
    (\code{\link{struc_sieve}}) are provided.}
  \item{core_args}{list of arguments for the core-generating function,
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
  \item{sub}{a character string used for plotting the subtitle.
    If \code{sub} is a logical and \code{TRUE} and \code{main} is unspecified, the
    name of the object supplied as \code{x} is used.}
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
  \item{title_margins}{either an object of class \code{"unit"} of length 2, or
    a numeric vector of length 2. The elements are recycled as needed.
    The two components specify the top and bottom \emph{title} margin
    of the plot, respectively. The default for each
    \emph{specified} title are 2 lines (and 0 else), except when a
    legend is plotted and \code{keep_aspect_ratio} is \code{TRUE}: in
    this case, the default values of both margins are set as to align
    the heights of legend and actual plot.
    When a numeric vector is supplied, the numbers are interpreted as
    \code{"lines"} units.  In addition, the unit or numeric vector
    may have named arguments (\samp{top} and \samp{bottom}), in which
    case the non-named argument specify the default value (recycled as
    needed), overloaded by the named arguments.}
  \item{legend_width}{An object of class \code{"unit"} of length
    1 specifying the width of the legend (if any). Default: 5 lines.}
  \item{pop}{logical indicating whether the generated viewport tree
    should be removed at the end of the drawing or not.}
  \item{main_gp, sub_gp}{object of class \code{"gpar"} containing the graphical
    parameters used for the main (sub) title, if specified.}
  \item{newpage}{logical indicating whether a new page should be created
    for the plot or not.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be
    fixed or not. If unspecified, the default is \code{TRUE} for
    two-dimensional tables and \code{FALSE} otherwise.}
  \item{prefix}{optional character string used as a prefix for the
    generated viewport and grob names.}
  \item{\dots}{For convenience, list of arguments passed to the
    labeling-generating function used.}
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

  Most elements of the plot, such as the core function, the spacing
  between the tiles, the shading of the tiles, the labeling, and the
  legend, are modularized in graphical appearance control (``grapcon'')
  functions and specified as parameters.  For
  each element \emph{foo} (= \code{spacing}, \code{labeling}, \code{core},
  or \code{legend}), \code{strucplot} takes two arguments:
  \var{foo} and \var{foo\_args}, which can be used to specify the
  parameters in the following alternative ways:

  \enumerate{
    \item Passing a suitable function to \var{foo} which subsequently
    will be called from \code{strucplot} to compute shadings, labelings,
    etc.
    \item Passing a corresponding \emph{generating} function to \var{foo},
    along with parameters passed to \var{foo\_args}, that generates such a
    function. Generating functions must inherit from classes
    \code{"grapcon_generator"} and \code{"}\var{foo}\code{"}.
    
    \item Except for the shading functions (\var{shading\_bar}),
    passing \var{foo(foo\_args)} to the \var{foo} argument.
    \item For shadings and spacings, passing the final parameter object
    itself; see the corresponding help pages for more details on the
    data structures.
  }
  
  If legends are drawn, a \sQuote{cinemascope}-like layout is used for
  the plot to preserve the 1:1 aspect ratio.

  If \code{type = "expected"}, the expected values are passed to the
  \code{observed} argument of the core function, and the observed
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

\references{
  Meyer D., Zeileis A., and Hornik K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.
}

\seealso{
  \code{\link{assoc}},
  \code{\link{mosaic}},
  \code{\link{sieve}},
  \code{\link{struc_assoc}},
  \code{\link{struc_sieve}},
  \code{\link{struc_mosaic}},
  \code{\link{structable}},
  \code{\link{doubledecker}},
  \code{\link{labelings}},
  \code{\link{shadings}},
  \code{\link{legends}},
  \code{\link{spacings}}
}
\examples{
data("Titanic")

strucplot(Titanic)
strucplot(Titanic, core = struc_assoc)
strucplot(Titanic, spacing = spacing_increase,
                   spacing_args = list(start = 0.5, rate = 1.5))
strucplot(Titanic, spacing = spacing_increase(start = 0.5, rate = 1.5))

## modify a tile's color
strucplot(Titanic, pop = FALSE)
grid.edit("rect:Class=1st,Sex=Male,Age=Adult,Survived=Yes",
          gp = gpar(fill = "red"))
}
\keyword{hplot}
