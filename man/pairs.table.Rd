\name{pairs.table}
\alias{pairs.table}
\title{Pairs Plot for Contingency Tables}
\description{
  Produces a matrix of strucplots.
}
\usage{
\method{pairs}{table}(x, upper_panel = pairs_mosaic, upper_panel_args = list(),
  lower_panel = pairs_mosaic, lower_panel_args = list(),
  diag_panel = pairs_barplot, diag_panel_args = list(),
  main = NULL, title_gp = gpar(fontsize = 20), space = 0.1,
  newpage = TRUE, \dots)
}
\arguments{
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute.}
  \item{upper_panel}{function for the upper triangle of the matrix, or
    corresponding generating function.}
  \item{upper_panel_args}{list of arguments for the generating
    function, if specified.}
  \item{lower_panel}{function for the lower triangle of the matrix, or
    corresponding generating function.}
  \item{lower_panel_args}{list of arguments for the panel-generating
    function, if specified.}
  \item{diag_panel}{function for the diagonal of the matrix, or
    corresponding generating function.}
  \item{diag_panel_args}{list of arguments for the generating
    function, if specified.}
  \item{main}{character string specifying the main title.}
  \item{title_gp}{object of class \code{"gpar"} used for the title.}
  \item{space}{double specifying the distance between the cells.}
  \item{newpage}{logical controlling whether a new grid page should be
    created.}
  \item{\dots}{currently not used.}
}
\details{
  This is a \code{\link[graphics]{pairs}} method for objects inheriting
  from class \code{"table"}.  It plots a matrix of pairwise mosaic
  plots.
  
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
  \code{\link{pairs_text}},
  \code{\link{pairs_barplot}},
  \code{\link{assoc}}
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
pairs(UCBAdmissions, diag_panel = pairs_text)
pairs(UCBAdmissions, upper_panel = pairs_assoc)
}
\keyword{hplot}
