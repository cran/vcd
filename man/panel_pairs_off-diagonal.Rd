\name{Pairs plot panel functions for off-diagonal cells}
\alias{pairs_strucplot}
\alias{pairs_mosaic}
\alias{pairs_assoc}
\title{Off-diagonal Panel Functions for Table Pairs Plot}
\description{
  Off-diagonal panel functions for \code{\link{pairs.table}}.
}
\usage{
pairs_strucplot(panel = mosaic, 
  type = c("pairwise", "total", "conditional", "joint"),
  legend = FALSE, margins = c(0, 0, 0, 0), labeling = NULL,
  shade = FALSE, \dots)
pairs_assoc(\dots)
pairs_mosaic(\dots)
}
\arguments{
  \item{panel}{function to be used for the plots in each
    cell---currently, \code{\link{assoc}} and \code{\link{mosaic}} are available.}
  \item{type}{character string specifying the type of independence model
    visualized in the cells.}
  \item{legend}{logical specifying whether a legend should be displayed
    in the cells or not.}
  \item{margins}{margins inside each cell (see \code{\link{strucplot}}).}
  \item{labeling}{labeling function or labeling-generating function (see
  \code{\link{strucplot}}).}
  \item{shade}{shading function or shading-generating function (see
  \code{\link{strucplot}}).}
  \item{\dots}{\code{pairs_mosaic} and \code{pairs_assoc}: parameters
    passed to \code{pairs_assocplot}. \code{pairs_assocplot}:
    other parameters passed to panel function.}
}
\details{
  These functions really just wrap \code{\link{assoc}} and
  \code{\link{mosaic}} by basically inhibiting labeling and
  legend-drawing and setting the margins to 0.
}
\value{
  A function with arguments:
  \item{x}{contingency table.}
  \item{i,j}{cell coordiantes.}
}
\seealso{
  \code{\link{pairs.table}},
  \code{\link{pairs_text}},
  \code{\link{pairs_barplot}},
  \code{\link{assoc}},
  \code{\link{mosaic}}
}
\references{
  Cohen, A. (1980),
  On the graphical display of the significant components in a two-way
  contingency table.
  \emph{Communications in Statistics---Theory and Methods}, \bold{A9},
  1025--1041.

  Friendly, M. (1992),
  Graphical methods for categorical data.
  \emph{SAS User Group International Conference Proceedings}, \bold{17},
  190--200.
  \url{http://www.math.yorku.ca/SCS/sugi/sugi17-paper.html}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data(UCBAdmissions)
data(PreSex)

pairs(PreSex)
pairs(UCBAdmissions)
pairs(UCBAdmissions, upper_panel_args = list(shade = FALSE))
pairs(UCBAdmissions, lower_panel = pairs_mosaic(type = "conditional"))
pairs(UCBAdmissions, upper_panel = pairs_assoc)
}
\keyword{hplot}
