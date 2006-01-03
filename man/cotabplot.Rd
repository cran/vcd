\name{cotabplot}
\alias{cotabplot}
\alias{cotabplot.default}
\alias{cotabplot.formula}
\title{Coplot for Contingency Tables}
\description{
  \code{cotabplot} is a generic function for creating trellis-like
  coplots (conditional plots) for contingency tables.
}
\usage{
cotabplot(x, \dots)
\method{cotabplot}{default}(x, cond = NULL,
  panel = cotab_mosaic, panel_args = list(),
  margins = rep(1, 4),
  text_gp = gpar(fontsize = 12), rect_gp = gpar(fill = grey(0.9)),
  pop = TRUE, newpage = TRUE,
  \dots)
\method{cotabplot}{formula}(formula, data = NULL, \dots)
}
\arguments{
  \item{x}{an object. The default method can deal with contingency tables
    in array form.}
  \item{cond}{margin index(es) or corresponding name(s) of the
    conditioning variables.}
  \item{panel}{panel function applied for each conditioned plot,
    see details.}
  \item{panel_args}{list of arguments passed to \code{panel} if this is
    a panel-generating function inheriting from class \code{"panel_generator"}.}
  \item{margins}{either an object of class \code{"unit"} of length 4, or
    a numeric vector of length 4. The elements are recycled as needed.
    giving the margins around the whole plot.}
  \item{text_gp}{object of class \code{"gpar"} used for the text in the panel titles.}
  \item{rect_gp}{object of class \code{"gpar"} used for the rectangles with the panel titles.}
  \item{pop}{logical indicating whether the generated viewport tree
    should be removed at the end of the drawing or not.}
  \item{newpage}{logical controlling whether a new grid page should be
    created.}
  \item{\dots}{further arguments passed to the panel-generating function.}
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}. It has to be of type \code{~ x + y | z}
    where \code{z} is/are the conditioning variable(s) used.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
}
\details{
  \code{cotabplot} is a generic function designed to create coplots or
  conditional plots (see Cleveland, 1993, and Becker, Cleveland, Shyu, 1996)
  similar to \code{\link{coplot}} but for contingency tables.
  
  \code{cotabplot} takes on computing the conditioning information 
  and setting up the trellis display, and then relies on a panel function
  to create plots from the full table and the conditioning information.
  A simple example would be a contingency table \code{tab} with margin
  names \code{"x"}, \code{"y"} and \code{"z"}. To produce this plot
  either the default interface can be used or the formula interface via
  
    \code{cotabplot(tab, "z")}
    \code{cotabplot(~ x + y | z, data = tab)}
    
  The panel function needs to be of the form
  
    \code{panel(x, condlevels)}
    
  where \code{x} is the \emph{full} table (\code{tab} in the example above)
  and \code{condlevels} is a named vector with the levels (e.g., 
  \code{c(z = "z1")} in the example above). 
  
  Alternatively, \code{panel} can also be a panel-generating function
  of class \code{"panel_generator"} which creates a function with the
  interface described above. The panel-generating function is called
  with the interface
  
    \code{panel(x, condvars, \dots)}
    
  where again \code{x} is the full table, \code{condvars} is now only
  a vector with the names of the conditioning variables (and not their
  levels, e.g., \code{"z"} in the example above). Further arguments
  can be passed to the panel-generating function via \code{\dots} which
  also includes the arguments set in \code{panel_args}.
  
  Suitable panel-generating functions for mosaic, association and sieve
  plots can be found at \code{\link{cotab_mosaic}}.
  
  A description of the underlying ideas is given in Zeileis, Meyer, Hornik
  (2005).
}
\seealso{
  \code{\link{cotab_mosaic}},
  \code{\link{cotab_coindep}},
  \code{\link{co_table}},
  \code{\link{coindep_test}}
}
\references{
  Becker, R.A., Cleveland, W.S., Shyu, M.-J. (1996),
  The visual design and control of trellis display.
  \emph{Journal of Computational and Graphical Statistics}, \bold{5},
  123--155.

  Cleveland, W.S. (1993), \emph{Visualizing Data}, Summit, New Jersey: Hobart Press.

  Zeileis, A., Meyer, D., Hornik K. (2005), \emph{Residual-based shadings for visualizing
  (conditional) independence}, Research Report, \bold{20}.
  \url{http://epub.wu-wien.ac.at/}
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
data("UCBAdmissions")

cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions)
cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions, panel = cotab_assoc)

ucb <- cotab_coindep(UCBAdmissions, condvars = "Dept", type = "assoc", n = 5000, margins = c(3, 1, 1, 3))
cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions, panel = ucb)
}
\keyword{hplot}
