\name{fourfold}
\alias{fourfold}
\title{Fourfold Plots}
\description{
  Creates an (extended) fourfold display of a \eqn{2 \times 2 \times k}{2 x 2 x k}
  contingency table, allowing for the visual inspection of the association
  between two dichotomous variables in one or several populations (strata).
}
\usage{
fourfold(x,
  color = c("#99CCFF", "#6699CC", "#FF5050", "#6060A0", "#FF0000", "#000080"),
  conf_level = 0.95, std = c("margins", "ind.max", "all.max"),
  margin = c(1, 2), space = 0.2, main = NULL, sub = NULL,
  mfrow = NULL, mfcol = NULL, extended = TRUE, ticks = 0.15,
  p_adjust_method = p.adjust.methods, newpage = TRUE,
  fontsize = 12)
}
\arguments{
  \item{x}{a \eqn{2 \times 2 \times k}{2 x 2 x k} contingency table in array form, or a
    \eqn{2 \times 2}{2 x 2} matrix if \eqn{k} is 1.}
  \item{color}{a vector of length 6 specifying the colors to use for the
    smaller and larger diagonals of each \eqn{2 \times 2}{2 x 2} table. The first pair
    is used for the standard (non-extended) plots, the other two
    for the extended version: the second/third pair is used for tables with
    non-significant/significant log-odds ratios, respectively, the latter
    being visualized in brighter colors.}
  \item{conf_level}{confidence level used for the confidence rings on
    the odds ratios.  Must be a single non-negative number less than 1;
    if set to 0, confidence rings are suppressed.}
  \item{std}{a character string specifying how to standardize the table.
    Must be one of \code{"margins"}, \code{"ind.max"}, or
    \code{"all.max"}, and can be abbreviated by the initial letter.
    If set to \code{"margins"}, each \eqn{2 \times 2}{2 x 2} table is standardized to
    equate the margins specified by \code{margin} while preserving the
    odds ratio.  If \code{"ind.max"} or \code{"all.max"}, the tables are
    either individually or simultaneously standardized to a maximal cell
    frequency of 1.}
  \item{margin}{a numeric vector with the margins to equate.  Must be
    one of \code{1}, \code{2}, or \code{c(1, 2)} (the default), which
    corresponds to standardizing only the row, only column, or both row and column in
    each \eqn{2 \times 2}{2 x 2} table.
    Only used if \code{std} equals \code{"margins"}.}
  \item{space}{the amount of space (as a fraction of the maximal radius
    of the quarter circles) used for the row and column labels.}
  \item{main, sub}{character string for the fourfold plot title/subtitle.}
  \item{mfrow, mfcol}{a numeric vector with two components:
    \var{nr} and \var{nc}, indicating that the displays for the \eqn{2 \times 2}{2 x 2}
    tables should be arranged in an \var{nr} by \var{nc} layout, filled
    by rows/columns.}
  \item{extended}{logical; if \code{TRUE}, extended plots are plotted,
    i.e., colors are brighter for significant log-odds ratios, and ticks
    are plotted for pies with positive log-odds.}
  \item{ticks}{the length of the ticks. If set to 0, no ticks are
    plotted.}
  \item{p_adjust_method}{method to be used for p-value adjustments for
    multi-stratum plots.  The p-values are used for the \sQuote{visual}
    significance tests of the odds ratios.}
  \item{newpage}{logical; if \code{TRUE}, \code{grid.newpage()} is
    called before plotting.}
  \item{fontsize}{fontsize of main title.}
}
\details{
  The fourfold display is designed for the display of \eqn{2 \times 2 \times k}{2 x 2 x k}
  tables.

  Following suitable standardization, the cell frequencies
  \eqn{f_{ij}}{f[i,j]} of each \eqn{2 \times 2}{2 x 2} table are shown as a quarter
  circle whose radius is proportional to
  \eqn{\sqrt{f_{ij}}}{sqrt(f[i,j])} so that its area is proportional to
  the cell frequency.  An association (odds ratio different from 1)
  between the binary row and column variables is indicated by the
  tendency of diagonally opposite cells in one direction to differ in
  size from those in the other direction; color is used to show this
  direction.  Confidence rings for the odds ratio allow a visual test of
  the null of no association; the rings for adjacent quadrants overlap
  iff the observed counts are consistent with the null hypothesis.

  Typically, the number \eqn{k} corresponds to the number of levels of a
  stratifying variable, and it is of interest to see whether the
  association is homogeneous across strata.  The fourfold display
  visualizes the pattern of association.  Note that the confidence rings
  for the individual odds ratios are not adjusted for multiple testing.
}
\references{
  Friendly, M. (1994),
  \emph{A fourfold display for 2 by 2 by \eqn{k} tables}.
  Technical Report 217, York University, Psychology Department,
  \url{http://www.math.yorku.ca/SCS/Papers/4fold/4fold.ps.gz}.

  Friendly, M. (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\seealso{
  \code{\link{mosaic}},
  \code{\link{assoc}}
}
\examples{
data("UCBAdmissions")
## Use the Berkeley admission data as in Friendly (1995).
x <- aperm(UCBAdmissions, c(2, 1, 3))
dimnames(x)[[2]] <- c("Yes", "No")
names(dimnames(x)) <- c("Sex", "Admit?", "Department")
ftable(x)

## Fourfold display of data aggregated over departments, with
## frequencies standardized to equate the margins for admission
## and sex.
## Figure 1 in Friendly (1994).
fourfold(margin.table(x, c(1, 2)))

## Fourfold display of x, with frequencies in each table
## standardized to equate the margins for admission and sex.
## Figure 2 in Friendly (1994).
fourfold(x)
cotabplot(x, panel = cotab_fourfold)

## Fourfold display of x, with frequencies in each table
## standardized to equate the margins for admission. but not
## for sex.
## Figure 3 in Friendly (1994).
fourfold(x, margin = 2)
}
\keyword{hplot}
