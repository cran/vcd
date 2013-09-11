\name{Pairs plot panel functions for diagonal cells}
\alias{pairs_barplot}
\alias{pairs_text}
\alias{pairs_diagonal_text}
\alias{pairs_diagonal_mosaic}
\title{Diagonal Panel Functions for Table Pairs Plot}
\description{
  Diagonal panel functions for \code{\link{pairs.table}}.
}
\usage{
pairs_barplot(gp_bars = NULL,
  gp_vartext = gpar(fontsize = 17),
  gp_leveltext = gpar(),
  just_leveltext = c("center", "bottom"),
  just_vartext = c("center", "top"),
  rot = 0, abbreviate = FALSE, check_overlap = TRUE, fill = "grey",
  var_offset = unit(1, "npc"), \dots)
pairs_text(dimnames = TRUE, gp_vartext = gpar(fontsize = 17), 
  gp_leveltext = gpar(), gp_border = gpar(), \dots) 
pairs_diagonal_text(varnames = TRUE, gp_vartext = gpar(fontsize = 17, fontface = "bold"),
  gp_leveltext = gpar(), gp_border = gpar(), pos = c("right","top"),
  distribute = c("equal","margin"), rot = 0, \dots)
pairs_diagonal_mosaic(split_vertical = TRUE, margins = unit(0, "lines"),
  offset_labels = -0.4, offset_varnames = 0,
  gp = NULL, fill = "grey", labeling = labeling_values, alternate_labels = TRUE, ...)
}
\arguments{
  \item{dimnames}{vector of logicals indicating whether the factor
    levels should be displayed (only used for \code{pairs_text}).}
  \item{varnames}{vector of logicals indicating whether the variable
    names should be displayed (only used for \code{pairs_text_diagonal}).}
  \item{gp_bars}{object of class \code{"gpar"} used for bars (only used
    for \code{pairs_barplot}). If unspecified, the default is to set the \code{fill}
    component of this object to the \code{fill} argument.}
  \item{gp_vartext}{object of class \code{"gpar"} used for the factor
    names.}
  \item{gp_leveltext}{object of class \code{"gpar"} used for the factor
    levels.}
  \item{gp_border}{object of class \code{"gpar"} used for the border
    (only used for \code{pairs_text}).}
  \item{gp}{object of class \code{"gpar"} used for the tiles (only used
    for \code{pairs_diagonal_mosaic}).
    If unspecified, the default is to set the \code{fill}
    component of this object to the \code{fill} argument.}
  \item{fill}{color vector or palette function used for the fill colors
    of bars (for \code{pairs_barplot}) or tiles (for
    \code{pairs_diagonal_mosaic}).}
  \item{labeling}{labeling function, passed to \code{mosaic()}}
  \item{alternate_labels}{should labels alternate top/bottom?}
  \item{just_leveltext, just_vartext}{character string indicating the
    justification for variable names and levels.}
  \item{pos}{character string of length 2 controlling the
    horizontal and vertical position of the variable names
    (only used for \code{pairs_text_diagonal}).}
  \item{rot}{rotation angle for the variable levels.}
  \item{distribute}{character string indicating whether levels should be
     distributed equally or according to the margins
     (only used for \code{pairs_text_diagonal}).}
  \item{abbreviate}{integer or logical indicating
    the number of characters the labels should be abbreviated
    to. \code{TRUE} means 1 character, \code{FALSE} causes no
    abbreviation.}
  \item{check_overlap}{If \code{TRUE}, some levels will suppressed to
    avoid overlapping, if any.}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where
    \eqn{k} is the number of margins of \code{x} (values are recycled as needed).
    A \code{TRUE} component indicates that the tile(s) of the
    corresponding dimension should be split vertically, \code{FALSE}
    means horizontal splits. Default is \code{FALSE.}}
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
  \item{offset_labels, offset_varnames}{numeric vector of length 4
    indicating the offset of the
    labels (variable names) for each of the four sides of the plot.}
  \item{var_offset}{object of class \code{"unit"} specifying the
        offset of variable names from the bottom of the bar plots
        created by \code{pairs_barplot}. If numeric, the unit defaults to "npc".}
  \item{\dots}{other parameters passed to the underlying graphics
    functions.}
}
\details{
  In the diagonal cells, the pairsplot visualizes statistics or
  information for each dimension (that is: the single factors) alone.
  \code{\link{pairs_text}} displays the factor's name, and optionally
  also the factor levels. \code{\link{pairs_barplot}} produces a bar plot
  of the corresponding factor, along with the factor's name.
}
\value{
  A function with one argument: the marginal table for the corresponding
  dimension.
}
\seealso{
  \code{\link{pairs.table}},
  \code{\link{pairs_assoc}},
  \code{\link{pairs_mosaic}}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data("UCBAdmissions")

pairs(UCBAdmissions) # pairs_barplot is default
pairs(UCBAdmissions, diag_panel = pairs_text)

pairs(UCBAdmissions, diag_panel = pairs_diagonal_text)
pairs(Titanic, diag_panel = pairs_diagonal_text)
pairs(Titanic, diag_panel = pairs_diagonal_text(distribute = "margin"))
pairs(Titanic,
      diag_panel = pairs_diagonal_text(distribute = "margin", rot = 45))
}
\keyword{hplot}
