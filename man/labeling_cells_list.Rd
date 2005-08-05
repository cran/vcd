\name{labeling_cells_list}
\alias{labeling_list}
\alias{labeling_cells}
\title{Labeling Functions for Strucplots}
\description{
  These functions generate labeling functions that produce labels for
  strucplots.
}
\usage{
labeling_cells(labels = TRUE, varnames = TRUE,
  abbreviate_labels = FALSE, abbreviate_varnames = FALSE,
  gp = gpar(), lsep = ": ", lcollapse = "\n",
  just = "center", pos = "center", rot = 0,
  margin = unit(0.5, "lines"), clip_cells = TRUE,
  text = NULL, \dots)
labeling_list(gp = gpar(), just = "left", pos = "left", lsep = ": ",
  sep = " ", offset = unit(c(2, 2), "lines"),
  varnames = TRUE, cols = 2, \dots) 
}
\arguments{
  \item{labels}{vector of logicals indicating, for each dimension,
    whether labels for the factor levels should be drawn or not.
    Values are recycled as needed.} 
  \item{varnames}{vector of logicals indicating, for each dimension,
    whether variable names should be drawn.
    Values are recycled as needed.}
  \item{abbreviate_labels}{vector of integers or logicals indicating,
    for each dimension, the number of characters the labels should be
    abbreviated to.
    \code{TRUE} means 1 character, \code{FALSE} causes no abbreviation.
    Values are recycled as needed.}
  \item{abbreviate_varnames}{vector of integers or logicals indicating,
    for each dimension, the number of characters the variable
    (i.e., dimension) names should be abbreviated to.
    \code{TRUE} means 1 character, \code{FALSE} causes no abbreviation.
    Values are recycled as needed.}
  \item{gp}{object of class \code{"gpar"} used for the text drawn.}
  \item{lsep}{character that separates variable names from the factor
    levels.}
  \item{sep}{character that separates the factor levels (only used for
    \code{labeling_list}).} 
  \item{offset}{object of class \code{"unit"} of length 2 specifying the
    offset in x- and y-direction of the text block drawn under the
    strucplot (only used for \code{labeling_list}).}
  \item{cols}{number of text columns (only used for
    \code{labeling_list}).}
  \item{lcollapse}{character that separates several variable name/factor
    level-combinations.  Typically a line break.
    (Only used for \code{labeling_cells}.)}
  \item{just, pos}{character string of length 1 (\code{labeling_list})
    or at most 2 (\code{labeling_cells}) specifying the labels'
    horizontal position and justification (horizontal and vertical for
    \code{labeling_cells}).}
  \item{rot}{rotation angle in degrees, used for all labels (only used
    for \code{labeling_cells}).}
  \item{margin}{object of class \code{"unit"} (a numeric value is
    converted to \code{"lines"}) specifying an offset from the cell
    borders (only used for \code{labeling_cells}).}
  \item{clip_cells}{logical indicating whether text should be clipped at
    the cell borders (only used for \code{labeling_cells}).}
  \item{text}{Optionally, a character table of the same dimensions than
    the contingency table whose entries will then be used instead of
    the labels.  \code{NA} entries are not drawn.
    This allows custom cell annotations (see examples).
    Only used for \code{labeling_cells}.}
  \item{\dots}{Currently not used.}
}
\details{
  These functions generate labeling functions that can add different
  kinds of labels to an existing plot. Typically they are 
  supplied to \code{\link{strucplot}} which then generates and calls
  the labeling function. They assume that a strucplot has been drawn
  and the corresponding viewport structure is pushed, so that by 
  navigating through the viewport tree the labels can be positioned
  appropriately.
  
  This help page only documents \code{labeling_list} and
  \code{labeling_cells}; more functions are described on the help page
  for \code{\link{labeling_border}}.
  
  The functions can also be used \sQuote{stand-alone} as shown in the
  examples.

  Using \code{labeling_list} will typically necessitate a bottom margin
  adjustment.
}
\value{
  A function with arguments:
  \item{d}{\code{"dimnames"} attribute from the visualized contingency
    table, or the visualized table itself from which the
    \code{"dimnames"} attributes will then be extracted.}
  \item{split_vertical}{vector of logicals indicating the split
    directions.}
  \item{condvars}{integer vector of conditioning dimensions}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{labeling_border}}, 
  \code{\link{structable}}, 
  \code{\link[grid]{grid.text}}
}
\examples{
mosaic(Titanic, labeling = labeling_cells)
mosaic(Titanic, labeling = labeling_list)

## A more complex example, adding the observed frequencies
## to a mosaic plot:
tab <- ifelse(Titanic < 6, NA, Titanic)
mosaic(Titanic, pop = FALSE)
labeling_cells(text = tab, margin = 0)(Titanic)
}
\keyword{hplot}

