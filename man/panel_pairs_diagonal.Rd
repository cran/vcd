\name{Pairs plot panel functions for diagonal cells}
\alias{pairs_barplot}
\alias{pairs_text}
\title{Diagonal Panel Functions for Table Pairs Plot}
\description{
  Diagonal panel functions for \code{\link{pairs.table}}.
}
\usage{
pairs_barplot(gp_bars = gpar(fill = "gray"),
  gp_vartext = gpar(fontsize = 17),
  gp_leveltext = gpar(), \dots)
pairs_text(dimnames = TRUE, gp_vartext = gpar(fontsize = 17), 
  gp_leveltext = gpar(), gp_border = gpar(), \dots) 
}
\arguments{
  \item{dimnames}{vector of logicals indicating whether the factor
    levels should be displayed (only used for \code{pairs_text}).}
  \item{gp_bars}{object of class \code{"gpar"} used for bars (only used
    for \code{pairs_barplot}).}
  \item{gp_vartext}{object of class \code{"gpar"} used for the factor
    names.}
  \item{gp_leveltext}{object of class \code{"gpar"} used for the factor
    levels.}
  \item{gp_border}{object of class \code{"gpar"} used for the border
    (only used for \code{pairs_text}).}
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
}
\keyword{hplot}
