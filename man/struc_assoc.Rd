\name{struc_assoc}
\alias{struc_assoc}
\title{Panel-generating Function for Association Plots}
\description{
  Panel-generating function for \code{strucplot} returning a function
  producing association plots.
}
\usage{
struc_assoc(compress = TRUE, xlim = NULL, ylim = NULL,
  yspace = unit(0.5, "lines"), xscale = 0.9, gp_axis = gpar(lty = 3))
}
\arguments{
  \item{compress}{logical; if \code{FALSE}, the space between the rows
    (columns) are chosen such that the \emph{total} heights (widths) of
    the rows (column) are all equal.  If \code{TRUE}, the space between
    the rows and columns is fixed and hence the plot is more
    \dQuote{compressed}.} 
  \item{xlim}{either a \eqn{2 \times k}{2 x k} matrix of doubles,
    \eqn{k} the number of total columns of the plot, or a recycled
    vector from which such a matrix will be constructed.
    The columns of \code{xlim} correspond to the columns of the
    association plot, the rows describe the column ranges (minimums in
    the first row, maximums in the second row).
    If \code{xlim} is \code{NULL}, the ranges are determined from the
    residuals according to \code{compress} (if \code{TRUE}: widest range
    from each column, if \code{FALSE}: from the whole association plot
    matrix).}
  \item{ylim}{either a \eqn{2 \times k}{2 x k} matrix of doubles,
    \eqn{k} the number of total rows of the plot, or a recycled vector
    from which such a matrix will be constructed.
    The columns of \code{ylim} correspond to the rows of the association
    plot, the rows describe the column ranges (minimums in the first
    row, maximums in the second row).
    If \code{ylim} is \code{NULL}, the ranges are determined from the
    residuals according to \code{compress} (if \code{TRUE}: widest range
    from each row, if \code{FALSE}: from the whole association plot
    matrix).}
  \item{xscale}{scale factor resizing the tile's width, thus adding
    additional space between the tiles. }
  \item{yspace}{object of class \code{"unit"} specifying additional
    space separating the rows.}
  \item{gp_axis}{object of class \code{"gpar"} specifying the visual
    aspects of the tiles' baseline.}
}
\details{
  This function is usually called by \code{strucplot} (typically when
  called by \code{assoc}) and returns a function used by
  \code{strucplot} to produce association plots. 
}
\value{
  A function with arguments:
  \item{residuals}{table of residuals.}
  \item{observed}{not used by \code{struc_assoc}.}
  \item{expected}{table of expected frequencies.}
  \item{spacing}{object of class \code{"unit"} specifying the space
    between the tiles.}
  \item{gp}{list of \code{gpar} objects used for the drawing the tiles.}
  \item{split_vertical}{vector of logicals indicating, for each
    dimension of the table, the split direction.}
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{assoc}},
  \code{\link{strucplot}},
  \code{\link{structable}}
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
\example{
## UCB Admissions
data(UCBAdmissions)
ucb <- aperm(UCBAdmissions)

## association plot for conditional independence
strucplot(ucb, expected = ~ Dept * (Admit + Gender),
  panel = struc_assoc(ylim = c(-4, 4)), labeling_args = list(abbreviate = c(Admit = 3)))
}
\keyword{hplot}
