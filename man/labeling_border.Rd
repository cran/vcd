\name{labeling_border}
\alias{labeling_border}
\alias{labeling_conditional}
\alias{labeling_left}
\alias{labeling_left2}
\alias{labeling_cboxed}
\alias{labeling_lboxed}
\alias{labeling_doubledecker}
\alias{labelings}
\title{Labeling Functions for Strucplots}
\description{
  These functions generate labeling functions used for strucplots.
}
\usage{
labeling_border(labels = TRUE, varnames = labels,
  set_labels = NULL, set_varnames = NULL, 
  tl_labels = NULL, alternate_labels = FALSE, tl_varnames = NULL, 
  gp_labels = gpar(fontsize = 12),
  gp_varnames = gpar(fontsize = 12, fontface = 2),
  rot_labels = c(0, 90, 0, 90), rot_varnames = c(0, 90, 0, 90),
  pos_labels = "center", pos_varnames = "center",
  just_labels = "center", just_varnames = pos_varnames,
  boxes = FALSE, fill_boxes = FALSE,
  offset_labels = c(0, 0, 0, 0), offset_varnames = offset_labels,
  labbl_varnames = NULL, labels_varnames = FALSE, sep = ": ",
  abbreviate = FALSE, rep = TRUE, clip = FALSE, \dots)
labeling_conditional(\dots)
labeling_left(rep = FALSE, pos_varnames = "left",
  pos_labels = "left", just_labels = "left", \dots)
labeling_left2(tl_labels = TRUE, clip = TRUE, pos_varnames = "left",
  pos_labels = "left", just_labels = "left", \dots)
labeling_cboxed(tl_labels = TRUE, boxes = TRUE, clip = TRUE,
  pos_labels = "center", \dots)
labeling_lboxed(tl_labels = FALSE, boxes = TRUE, clip = TRUE,
  pos_labels = "left", just_labels = "left",
  labbl_varnames = FALSE, \dots)
labeling_doubledecker(lab_pos = c("bottom", "top"), dep_varname = TRUE,
  boxes = NULL, clip = NULL, labbl_varnames = FALSE,
  rot_labels = rep.int(0, 4),
  pos_labels = c("left", "center", "left", "center"),
  just_labels = c("left", "left", "left", "center"),
  varnames = NULL, gp_varnames = gpar(fontsize = 12, fontface = 2),
  offset_varnames = c(0, -0.6, 0, 0), tl_labels = NULL, \dots)
}
\arguments{
  \item{labels}{vector of logicals indicating whether labels should be
    drawn for a particular dimension.}
  \item{varnames}{vector of logicals indicating whether variable names
    should be drawn for a particular dimension.}
  \item{set_labels}{An optional character vector with named components
    replacing the so-specified variable names. The component names must
    exacty match the variable names to be replaced.}
  \item{set_varnames}{An optional list with named components of
    character vectors replacing the labels of the so-specified
    variables. The component names must
    exacty match the variable names whose labels should be replaced.}
  \item{tl_labels}{vector of logicals indicating whether labels should
    be positioned on top (column labels) / left (row labels) for a
    particular dimension.}
  \item{alternate_labels}{vector of logicals indicating whether labels should
    be alternated on the top/bottom (left/right) side of the plot for a
    particular dimension.}
  \item{tl_varnames}{vector of logicals indicating whether variable
    names should be positioned on top (column labels) / on left (row
    labels) for a particular dimension.}
  \item{gp_labels}{list of objects of class \code{"gpar"} used for
    drawing the labels.}
  \item{gp_varnames}{list of objects of class \code{"gpar"} used for
    drawing the variable names.}
  \item{rot_labels}{vector of rotation angles for the labels for each of
    the four sides of the plot.}
  \item{rot_varnames}{vector of rotation angles for the variable names
    for each of the four sides of the plot.}
  \item{pos_labels}{character string of label positions (\code{"left"},
    \code{"center"}, \code{"right"}) for each of the variables.}
  \item{pos_varnames}{character string of variable names positions
    (\code{"left"}, \code{"center"}, \code{"right"}) for each of the
    four sides of the plot.}
  \item{just_labels}{character string of label justifications
    (\code{"left"}, \code{"center"}, \code{"right"}) for each of the
    variables.}
  \item{just_varnames}{character string of variable names justifications
    (\code{"left"}, \code{"center"}, \code{"right"}) for each of the
    four sides of the plot.}
  \item{boxes}{vector of logicals indicating whether boxes should be
    drawn around the labels for a particular dimension.}
  \item{fill_boxes}{Either a vector of logicals, or a vector of characters,
    or a list of such vectors, specifying the fill colors for the
    boxes. \code{"TRUE"} and \code{"FALSE"} values are
    transformed into \code{"grey"} and \code{"white"}, respectively. If
    \code{fill_boxes} is atomic, each component
    specifies a basic color for the corresponding dimension. This color
    is transformed into its HSV representation, and the value is varied
    from 50\% to 100\% to give a sequential color palette for the
    levels. For \code{NA} components, no palette is produced (no fill
    color). If \code{fill_boxes} is a list of vectors,
    each vector specifies the level colors of the corresponding
    dimension.}
  \item{offset_labels, offset_varnames}{numeric vector of length 4
    indicating the offset of the
    labels (variable names) for each of the four sides of the plot.}
  \item{labbl_varnames}{vector of logicals indicating whether variable
    names should be drawn on the left (column variables) / on top (row
    variables) of the corresponding labels.}
  \item{labels_varnames}{vector of logicals indicating, for each
    dimension, whether the variable name should be added to the
    corresponding labels or not.}
  \item{sep}{separator used if any component of \code{"labels_varnames"}
    is \code{TRUE}.}
  \item{abbreviate}{vector of integers or logicals indicating, for each
    dimension, the number of characters the labels should be abbreviated
    to.
    \code{TRUE} means 1 character, \code{FALSE} causes no abbreviation.
    Values are recycled as needed.}
  \item{rep}{vector of logicals indicating, for each dimension, whether
    labels should be repeated for all conditioning strata, or appear
    only once.}
  \item{clip}{vector of integers indicating, for each dimension, whether
    labels should be clipped to not overlap.}
  \item{lab_pos}{character string switching between \code{"top"} or
    \code{"bottom"} position of the labels (only used for
    \code{labeling_doubledecker}).}
  \item{dep_varname}{logical or character string. If logical, this is
    indicating whether the name of the dependent varable should be
    printed or not. A character string will be printed instead of the
    variable name taken from the dimnames.}
  \item{\dots}{only used for \code{labeling_conditional} and
    \code{labeling_doubledecker}: parameters
    passed to \code{labeling_cells} and \code{labeling_border}.}
}
\details{
  These functions generate labeling functions called by
  \code{\link{strucplot}} for their side-effect of adding labels to the
  plot. They suppose that a strucplot has been drawn and the
  corresponding viewport structure is pushed, since the positions of the
  viewports are used for the label positioning.
  Note that the functions can also be used \sQuote{stand-alone} as
  shown in the examples.
  
  All values supplied to vectorized arguments can be
  \sQuote{abbreviated} by using named components which override the
  default component values.  In addition, these defaults can be
  overloaded by the sequence of non-named components which are recycled
  as needed (see examples).

  This help page only documents \code{labeling_border} and
  derived functions, more functions are described on the help page
  for \code{\link{labeling_cells}} and \code{\link{labeling_list}}.

  \code{labeling_left}, \code{labeling_left2}, \code{labeling_cboxed},
  and \code{labeling_lboxed}
  are really just wrappers to \code{labeling_border}, and good examples for
  the parameter usage.
}
\value{
  A function with arguments:
  \item{d}{\code{"dimnames"} attribute from the visualized contingency
    table, or the visualized table itself from which the
    \code{"dimnames"} attributes will then be extracted.}
  \item{split_vertical}{vector of logicals indicating the split
    directions.}
  \item{condvars}{integer vector of conditioning dimensions.}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\references{
  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.
}
  
\seealso{
  \code{\link{labeling_cells}}, 
  \code{\link{labeling_list}}, 
  \code{\link{structable}}, 
  \code{\link[grid]{grid.text}}
}
\examples{
data("Titanic")

mosaic(Titanic)
mosaic(Titanic, labeling = labeling_left)
labeling_left

mosaic(Titanic, labeling = labeling_cboxed)
labeling_cboxed

mosaic(Titanic, labeling = labeling_lboxed)
labeling_lboxed

data("PreSex")
mosaic(~ PremaritalSex + ExtramaritalSex | Gender + MaritalStatus,
  data = PreSex, labeling = labeling_conditional)

## specification of vectorized arguments 
mosaic(Titanic, labeling_args = list(abbreviate = c(Survived = TRUE)))

mosaic(Titanic, labeling_args = list(clip = TRUE, boxes = TRUE,
  fill_boxes = c(Survived = "green", "red")))

mosaic(Titanic, labeling_args = list(clip = TRUE, boxes = TRUE,
  fill_boxes = list(Sex = "red", "green")))

mosaic(Titanic, labeling_args = list(clip = TRUE, boxes = TRUE,
  fill_boxes = list(Sex = c(Male = "red", "blue"), "green")))

## change variable names
mosaic(Titanic, labeling_args = list(set_varnames = c(Sex = "Gender")))

## change labels
mosaic(Titanic, labeling_args = list(set_varnames = c(Survived = "Status"),
  set_labels = list(Survived = c("Survived", "Not Survived")), rep = FALSE))
}
\keyword{hplot}

