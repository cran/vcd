\name{cotab_panel}
\alias{cotab_mosaic}
\alias{cotab_assoc}
\alias{cotab_sieve}
\alias{cotab_coindep}
\alias{cotab_fourfold}
\title{Panel-generating Functions for Contingency Table Coplots}
\description{
  Panel-generating functions visualizing contingency tables that
  can be passed to \code{cotabplot}.
}
\usage{
cotab_mosaic(x = NULL, condvars = NULL, \dots)
cotab_assoc(x = NULL, condvars = NULL, ylim = NULL, \dots)
cotab_sieve(x = NULL, condvars = NULL, \dots)
cotab_fourfold(x = NULL, condvars = NULL, \dots)
cotab_coindep(x, condvars,
  test = c("doublemax", "maxchisq", "sumchisq"),
  level = NULL, n = 1000, interpolate = c(2, 4),
  h = NULL, c = NULL, l = NULL, lty = 1,
  type = c("mosaic", "assoc"), legend = FALSE, ylim = NULL, \dots)
}
\arguments{
  \item{x}{a contingency tables in array form.}
  \item{condvars}{margin name(s) of the conditioning variables.}
  \item{ylim}{y-axis limits for \code{assoc} plot. By default this
    is computed from \code{x}.}
  \item{test}{character indicating which type of statistic
    should be used for assessing conditional independence.}
  \item{level,n,h,c,l,lty,interpolate}{variables controlling the HCL shading of the
    residuals, see \code{\link{shadings}} for more details.}
  \item{type}{character indicating which type of plot should be produced.}
  \item{legend}{logical. Should a legend be produced in each panel?}
  \item{\dots}{further arguments passed to the plotting function (such as
    \code{\link{mosaic}} or \code{\link{assoc}} or \code{\link{sieve}}
    respectively).}
}
\details{
  These functions of class \code{"panel_generator"} are panel-generating
  functions for use with \code{\link{cotabplot}}, i.e., they return functions
  with the interface
  
    \code{panel(x, condlevels)}

  required for \code{cotabplot}. The functions produced by \code{cotab_mosaic},
  \code{cotab_assoc} and \code{cotab_sieve} essentially only call \code{co_table}
  to produce the conditioned table and then call \code{\link{mosaic}}, \code{\link{assoc}}
  or \code{\link{sieve}} respectively with the arguments specified.
  
  The function \code{cotab_coindep} is similar but additionally chooses an appropriate
  residual-based shading visualizing the associated conditional independence
  model. The conditional independence test is carried out via \code{\link{coindep_test}}
  and the shading is set up via \code{\link{shading_hcl}}.
    
  A description of the underlying ideas is given in Zeileis, Meyer, Hornik
  (2005).
}
\seealso{
  \code{\link{cotabplot}},
  \code{\link{mosaic}},
  \code{\link{assoc}},
  \code{\link{sieve}},
  \code{\link{co_table}},
  \code{\link{coindep_test}},
  \code{\link{shading_hcl}}
}
\references{
  Meyer, D., Zeileis, A., and Hornik, K. (2006),
  The strucplot framework: Visualizing multi-way contingency tables with
  \pkg{vcd}.
  \emph{Journal of Statistical Software}, \bold{17(3)}, 1-48.
  URL http://www.jstatsoft.org/v17/i03/ and available as 
  \code{vignette("strucplot")}.

  Zeileis, A., Meyer, D., Hornik K. (2007), \emph{Residual-based shadings for visualizing
  (conditional) independence}, \emph{Journal of Computational and Graphical Statistics}, \bold{16},
  507--525.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\examples{
data("UCBAdmissions")

cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions)
cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions, panel = cotab_assoc)
cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions, panel = cotab_fourfold)

ucb <- cotab_coindep(UCBAdmissions, condvars = "Dept", type = "assoc", n = 5000, margins = c(3, 1, 1, 3))
cotabplot(~ Admit + Gender | Dept, data = UCBAdmissions, panel = ucb)
}
\keyword{hplot}
