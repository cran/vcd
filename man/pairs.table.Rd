\name{pairs.table}
\alias{pairs.table}
\alias{pairs.structable}
\title{Pairs Plot for Contingency Tables}
\description{
  Produces a matrix of strucplot displays.
}
\usage{
\method{pairs}{table}(x, upper_panel = pairs_mosaic, upper_panel_args = list(),
  lower_panel = pairs_mosaic, lower_panel_args = list(),
  diag_panel = pairs_barplot, diag_panel_args = list(),
  main = NULL, sub = NULL, main_gp = gpar(fontsize = 20),
  sub_gp = gpar(fontsize = 15), space = 0.3,
  newpage = TRUE, pop = TRUE, margins = unit(1, "lines"), \dots)
}
\arguments{

  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute.}
  \item{upper_panel}{function for the upper triangle of the matrix, or
    corresponding generating function. If \code{NULL}, no panel is drawn.}
  \item{upper_panel_args}{list of arguments for the generating
    function, if specified.}
  \item{lower_panel}{function for the lower triangle of the matrix, or
    corresponding generating function. If \code{NULL}, no panel is drawn.}
  \item{lower_panel_args}{list of arguments for the panel-generating
    function, if specified.}
  \item{diag_panel}{function for the diagonal of the matrix, or
    corresponding generating function. If \code{NULL}, no panel is drawn.}
  \item{diag_panel_args}{list of arguments for the generating
    function, if specified.}
  \item{main}{either a logical, or a character string used for plotting
    the main title.  If \code{main} is a logical and \code{TRUE}, the
    name of the object supplied as \code{x} is used.}
  \item{sub}{a character string used for plotting the subtitle.
    If \code{sub} is a logical and \code{TRUE} and \code{main} is unspecified, the
    name of the object supplied as \code{x} is used.}
  \item{main_gp, sub_gp}{object of class \code{"gpar"} containing the graphical
    parameters used for the main (sub) title, if specified.}
  \item{space}{double specifying the distance between the cells.}
  \item{newpage}{logical controlling whether a new grid page should be
    created.}
  \item{pop}{logical indicating whether all viewports should be popped
    after the plot has been drawn.}
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
  \item{\dots}{For convenience,
    list of arguments for the panel-generating
    functions of upper and lower panels, if specified.}
}
\details{
  This is a \code{\link[graphics]{pairs}} method for objects inheriting
  from class \code{"table"} or \code{"structable"}.
  It plots a matrix of pairwise mosaic plots.
  
  Four independence types are distinguished: \code{"pairwise"},
  \code{"total"}, \code{"conditional"} and \code{"joint"}.
  The pairwise mosaic matrix shows bivariate marginal relations,
  collapsed over all other variables.
  The total independence mosaic matrix shows mosaic plots for mutual
  independence, i.e., for marginal and conditional independence among
  all pairs of variables.
  The conditional independence mosaic matrix shows mosaic plots for
  marginal independence given all other variables.
  The joint independence mosaic matrix shows mosaic plots for joint
  independence of all pairs of variables from the others.

  This method uses panel functions called for each cell of the
  matrix which can be different for upper matrix, lower matrix, and
  diagonal cells.  Correspondingly, for each panel parameter \var{foo}
  (= \samp{upper}, \samp{lower}, or \samp{diag}), \code{pairs.table} takes
  two arguments: \var{foo\_panel} and \var{foo\_panel\_args}, which can
  be used to specify the parameters as follows:

  \enumerate{
    \item Passing a suitable panel function to \var{foo\_panel} which
    subsequently is called for each cell with the corresponding
    coordinates.
    \item Passing a corresponding \emph{generating function} (of class
    \code{"panel_generator"}) to \var{foo\_panel}, along with parameters passed to
    \var{foo\_panel\_args}, that generates such a function. 
  }
Hence, the second approach is equivalent to the first if
\var{foo\_panel(foo\_panel\_args)} is passed to
\var{foo\_panel}.

  
}
\seealso{
  \code{\link{pairs_mosaic}},
  \code{\link{pairs_assoc}},
  \code{\link{pairs_sieve}},
  \code{\link{pairs_diagonal_text}},
  \code{\link{pairs_text}},
  \code{\link{pairs_barplot}},
  \code{\link{assoc}},
  \code{\link{sieve}},
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
data("UCBAdmissions")
data("PreSex")
data(HairEyeColor)
hec = structable(Eye ~ Sex + Hair, data = HairEyeColor)

pairs(PreSex)
pairs(UCBAdmissions)
pairs(UCBAdmissions, upper_panel_args = list(shade = TRUE))
pairs(UCBAdmissions, lower_panel = pairs_mosaic(type = "conditional"))
pairs(UCBAdmissions, diag_panel = pairs_text)
pairs(UCBAdmissions, upper_panel = pairs_assoc, shade = TRUE)
pairs(hec, highlighting = 2, diag_panel_args = list(fill = grey.colors))
pairs(hec, highlighting = 2, diag_panel = pairs_diagonal_mosaic,
           diag_panel_args = list(fill = grey.colors, alternate_labels =TRUE))

}
\keyword{hplot}
