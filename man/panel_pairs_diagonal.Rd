\name{Pairs plot panel functions for diagonal cells}
\alias{pairs_barplot}
\alias{pairs_text}
\alias{pairs_diagonal_text}
\title{Diagonal Panel Functions for Table Pairs Plot}
\description{
  Diagonal panel functions for \code{\link{pairs.table}}.
}
\usage{
pairs_barplot(gp_bars = gpar(fill = "gray"),
  gp_vartext = gpar(fontsize = 17),
  gp_leveltext = gpar(),
  just_leveltext = c("center", "bottom"),
  just_vartext = c("center", "top"),
  rot = 0, abbreviate = FALSE, \dots)
pairs_text(dimnames = TRUE, gp_vartext = gpar(fontsize = 17), 
  gp_leveltext = gpar(), gp_border = gpar(), \dots) 
pairs_diagonal_text(varnames = TRUE, gp_vartext = gpar(fontsize = 17, fontface = "bold"),
  gp_leveltext = gpar(), gp_border = gpar(), pos = c("right","top"),
  distribute = c("equal","margin"), rot = 0, \dots) 
}
\arguments{
  \item{dimnames}{vector of logicals indicating whether the factor
    levels should be displayed (only used for \code{pairs_text}).}
  \item{varnames}{vector of logicals indicating whether the variable
    names should be displayed (only used for \code{pairs_text_diagonal}).}
  \item{gp_bars}{object of class \code{"gpar"} used for bars (only used
    for \code{pairs_barplot}).}
  \item{gp_vartext}{object of class \code{"gpar"} used for the factor
    names.}
  \item{gp_leveltext}{object of class \code{"gpar"} used for the factor
    levels.}
  \item{gp_border}{object of class \code{"gpar"} used for the border
    (only used for \code{pairs_text}).}
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
    to. \code{TRUE} means 1 character, \code{FALSE} causes no abbreviation.}
  \item{\dots}{other parameters passed to the underlying graphics
    functions.}
}
\details{
  In the diagonal cells, the pairsplot visualizes statistics or
  information for each dimension (that is: the single factors) alone.
  \code{\link{pairs_text}} displays the factor's name, and optionally
  also the factor levels. \code{\link{pairs_barplot}} produces a barplot
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
data(UCBAdmissions)

pairs(UCBAdmissions) # pairs_barplot is default
pairs(UCBAdmissions, diag_panel = pairs_text)

pairs(UCBAdmissions, diag_panel = pairs_diagonal_text)
pairs(Titanic, diag_panel = pairs_diagonal_text)
pairs(Titanic, diag_panel = pairs_diagonal_text(distribute = "margin"))
pairs(Titanic,
      diag_panel = pairs_diagonal_text(distribute = "margin", rot = 45))
}
\keyword{hplot}
