\name{assoc}
\alias{assoc}
\alias{assoc.default}
\alias{assoc.formula}
\title{Extended Association Plots}
\description{
  Produce an association plot indicating deviations from a
  specified independence model in a possibly high-dimensional
  contingency table.
}
\usage{
\method{assoc}{default}(x, row_vars = NULL, col_vars = NULL, compress = TRUE,
  xlim = NULL, ylim = NULL,
  spacing = spacing_conditional(sp = 0), spacing_args = list(),
  split_vertical = NULL, keep_aspect_ratio = FALSE, 
  xscale = 0.9, yspace = unit(0.5, "lines"), main = NULL, sub = NULL,
  \dots, residuals_type = "Pearson", gp_axis = gpar(lty = 3))
\method{assoc}{formula}(formula, data = NULL, \dots, subset = NULL, na.action = NULL, main = NULL, sub = NULL)
}
\arguments{
  \item{x}{a contingency table in array form with optional category
    labels specified in the \code{dimnames(x)} attribute, or an object
    inheriting from the \code{"ftable"} class (such as
    \code{"structable"} objects).}
  \item{row_vars}{a vector of integers giving the indices, or a
    character vector giving the names of the variables to be used for
    the rows of the association plot.}
  \item{col_vars}{a vector of integers giving the indices, or a
    character vector giving the names of the variables to be used for
    the columns of the association plot.}
  \item{compress}{logical; if \code{FALSE}, the space between the rows
    (columns) are chosen such that the \emph{total} heights (widths) of
    the rows (columns) are all equal.  If \code{TRUE}, the space between
    rows and columns is fixed and hence the plot is more
    \dQuote{compressed}.}
  \item{xlim}{a \eqn{2 \times k}{2 x k} matrix of doubles, \eqn{k}
    number of total columns of the plot.  The columns of \code{xlim}
    correspond to the columns of the association plot, the rows describe
    the column ranges (minimums in the first row, maximums in the second
    row).  If \code{xlim} is \code{NULL}, the ranges are determined from
    the residuals according to \code{compress} (if \code{TRUE}: widest
    range from each column, if \code{FALSE}: from the whole association
    plot matrix).}
  \item{ylim}{a \eqn{2 \times k}{2 x k} matrix of doubles, \eqn{k}
    number of total rows of the plot.  The columns of \code{ylim}
    correspond to the rows of the association plot, the rows describe
    the column ranges (minimums in the first row, maximums in the second
    row).  If \code{ylim} is \code{NULL}, the ranges are determined from
    the residuals according to \code{compress} (if \code{TRUE}: widest
    range from each row, if \code{FALSE}: from the whole association
    plot matrix).}
  \item{spacing}{a spacing object, a spacing function, or a
    corresponding generating function (see \code{\link{strucplot}} for
    more information).  The default is the spacing-generating function
    \code{\link{spacing_conditional}} that is (by default) called with the
    argument list \code{spacing_args} (see \code{spacings} for more details).}
  \item{spacing_args}{list of arguments for the spacing-generating function, if
    specified (see \code{\link{strucplot}} for more information).}
  \item{split_vertical}{vector of logicals of length \eqn{k}, where \eqn{k}
    is the number of margins of \code{x} (default: \code{FALSE}).
    Values are recycled as needed.
    A \code{TRUE} component indicates that the corresponding dimension
    is folded into the columns, \code{FALSE} folds the dimension into the
    rows.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be
    fixed or not.}
  \item{residuals_type}{a character string indicating the type of
    residuals to be computed. Currently, only Pearson residuals are
    supported.}
  \item{xscale}{scale factor resizing the tile's width, thus adding
    additional space between the tiles. }
  \item{yspace}{object of class \code{"unit"} specifying additional
    space separating the rows.}
  \item{gp_axis}{object of class \code{"gpar"} specifying the visual
    aspects of the tiles' baseline.}
  \item{formula}{a formula object with possibly both left and right hand
    sides specifying the column and row variables of the flat table.}
  \item{data}{a data frame, list or environment containing the variables
    to be cross-tabulated, or an object inheriting from class \code{table}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used. Ignored if \code{data} is a contingency table.}
  \item{na.action}{an optional function which indicates what should happen when
    the data contain \code{NA}s. Ignored if \code{data} is a contingency table.}
  \item{main, sub}{either a logical, or a character string used for plotting
    the main (sub) title.  If logical and \code{TRUE}, the
    name of the \code{data} object is used.}
  \item{\dots}{other parameters passed to \code{\link{strucplot}}}
}
\details{
  Association plots have been suggested by Cohen (1980) and extended
  by Friendly (1992) and provide a means for visualizing the residuals
  of an independence model for a contingency table.

  \code{assoc} is a generic function and currently has a default method and a
  formula interface.  Both are high-level interfaces to the
  \code{\link{strucplot}} function, and produce (extended) association
  plots.  Most of the functionality is described there, such as
  specification of the independence model, labeling, legend, spacing,
  shading, and other graphical parameters.

  For a contingency table, the signed contribution
  to Pearson's \eqn{\chi^2}{chi^2} for cell \eqn{\{ij\ldots k\}} is
  
  \deqn{d_{ij\ldots k} = \frac{(f_{ij\ldots k} - e_{ij\ldots k})}{
      \sqrt{e_{ij\ldots k}}}}{d_\{ij\ldotsk\} = (f_\{ij\ldotsk\} -
    e_\{ij\ldotsk\}) / sqrt(e_\{ij\ldotsk\})}
  
  where \eqn{f_{ij\ldots k}}{f_\{ij\ldotsk\}} and \eqn{e_{ij\ldots
      k}}{e_\{ij\ldotsk\}} 
  are the observed and  expected counts corresponding to the cell.  In
  the association plot, each cell is represented by a
  rectangle that has (signed) height proportional to \eqn{d_{ij\ldots
      k}}{d_\{ij\ldotsk\}} 
  and width proportional to
  \eqn{\sqrt{e_{ij\ldots k}}}{sqrt(e_\{ij...k\})},
  so that the area of the box is proportional to the difference in
  observed and expected frequencies.  The rectangles in each row are
  positioned relative to a baseline indicating independence
  (\eqn{d_{ij\ldots k} = 0}{d_\{ij\ldotsk\} = 0}). 
  If the observed frequency of a cell is greater than the expected one,
  the box rises above the baseline, and falls below otherwise.

  Additionally, the residuals can be colored depending on a specified
  shading scheme (see Meyer et al., 2003).  Package \pkg{vcd} offers a range of
  \emph{residual-based} shadings (see the shadings help page). Some of
  them allow, e.g., the visualization of test statistics.

  Unlike the \code{\link[graphics]{assocplot}} function in the
  \pkg{graphics} package, this function allows the visualization of
  contingency tables with more than two dimensions. Similar to the
  construction of \sQuote{flat} tables (like objects of class \code{"ftable"} or
  \code{"structable"}), the dimensions are folded into rows and columns.
 
  The layout is very flexible: the specification of shading, labeling,
  spacing, and legend is modularized (see \code{\link{strucplot}} for
  details).

}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\seealso{
  \code{\link{mosaic}},
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
  \url{http://datavis.ca/papers/sugi/sugi17.pdf}
  
  Meyer, D., Zeileis, A., Hornik, K. (2003), Visualizing independence using
  extended association plots. \emph{Proceedings of the 3rd International
  Workshop on Distributed Statistical Computing}, K. Hornik, F. Leisch,
  A. Zeileis (eds.), ISSN 1609-395X.
  \url{https://www.R-project.org/conferences/DSC-2003/Proceedings/}

  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.

}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\examples{
data("HairEyeColor")
## Aggregate over sex:
(x <- margin.table(HairEyeColor, c(1, 2)))

## Ordinary assocplot:
assoc(x)
## and with residual-based shading (of independence)
assoc(x, main = "Relation between hair and eye color", shade = TRUE)

## Aggregate over Eye color:
(x <- margin.table(HairEyeColor, c(1, 3)))
chisq.test(x)
assoc(x, main = "Relation between hair color and sex", shade = TRUE)

# Visualize multi-way table
assoc(aperm(HairEyeColor), expected = ~ (Hair + Eye) * Sex,
      labeling_args = list(just_labels = c(Eye = "left"),
                           offset_labels = c(right = -0.5),
                           offset_varnames = c(right = 1.2),
                           rot_labels = c(right = 0),
                           tl_varnames = c(Eye = TRUE))
)

assoc(aperm(UCBAdmissions), expected = ~ (Admit + Gender) * Dept, compress = FALSE,
      labeling_args = list(abbreviate = c(Gender = TRUE), rot_labels = 0)
)
}
\keyword{hplot}
