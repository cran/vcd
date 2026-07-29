\name{rootogram}
\alias{rootogram}
\alias{rootogram.default}
\alias{rootogram.goodfit}
\title{Rootograms}
\description{
  Rootograms of observed and fitted values.
}
\usage{
\method{rootogram}{default}(x, fitted, names = NULL, scale = c("sqrt", "raw"),
  type = c("hanging", "standing", "deviation"),
  shade = FALSE, legend = TRUE,
  legend_args = list(x = 0, y = 0.2, height = 0.6), df = NULL,
  rect_gp = NULL, rect_gp_args = list(),
  lines_gp = gpar(col = "red", lwd = 2),
  points_gp = gpar(col = "red"), pch = 19,
  xlab = NULL, ylab = NULL, ylim = NULL,
  main = NULL, sub = NULL,
  margins = unit(0, "lines"),
  title_margins = NULL, legend_width = NULL,
  main_gp = gpar(fontsize = 20),
  sub_gp = gpar(fontsize = 15),
  name = "rootogram", prefix = "",
  keep_aspect_ratio = FALSE,
  newpage = TRUE, pop = TRUE,
  return_grob = FALSE, \dots)
}
\arguments{
  \item{x}{either a vector or a 1-way table of frequencies.}
  \item{fitted}{a vector of fitted frequencies.}
  \item{names}{a vector of names passed to \code{\link{grid_barplot}}, if set
    to \code{NULL} the names of \code{x} are used.}
  \item{scale}{a character string indicating whether the values should be
    plotted on the raw or square root scale.}
  \item{type}{a character string indicating if the bars for the observed
    frequencies should be \code{hanging} or \code{standing} or indicate the \code{deviation}
    between observed and fitted frequencies.}

  \item{shade}{logical specifying whether \code{rect_gp} should be set
    to colors corresponding to the pearson residuals, i.e., if a
    residual-based shading should be applied to the bars.}
  
  \item{legend}{either a legend-generating function, or a legend
    function (see details and \code{\link{legends}}), or a logical.
    If \code{legend} is \code{NULL} or \code{TRUE} and \code{gp} is a
    function, legend defaults to \code{\link{legend_resbased}}.}
  \item{legend_args}{list of arguments for the legend-generating
    function, if specified.}
  \item{df}{degrees of freedom passed to the shading functions
    used for inference.}
  
  \item{rect_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the rectangles, shading function or a
    corresponding generating function (see \code{\link{shadings}}). If
    unspecified and no shading is applied, defaults to light grey fill color for the bars.}
  \item{rect_gp_args}{list of arguments for the shading-generating
    function, if specified for \code{rect_gp}.}
  \item{lines_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the lines.}
  \item{points_gp}{a \code{"gpar"} object controlling the grid graphical
    parameters of the points.}
  \item{pch}{plotting character for the points.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{ylim}{limits for the y axis.}
  \item{main}{either a logical, or a character string used for plotting
    the main title.  If \code{main} is a logical and \code{TRUE}, the
    name of the object supplied as \code{x} is used.}
  \item{sub}{a character string used for plotting the subtitle.
    If \code{sub} is a logical and \code{TRUE} and \code{main} is unspecified, the
    name of the object supplied as \code{x} is used.}
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
  \item{title_margins}{either an object of class \code{"unit"} of length 2, or
    a numeric vector of length 2. The elements are recycled as needed.
    The two components specify the top and bottom \emph{title} margin
    of the plot, respectively. The default for each
    \emph{specified} title are 2 lines (and 0 else), except when a
    legend is plotted and \code{keep_aspect_ratio} is \code{TRUE}: in
    this case, the default values of both margins are set as to align
    the heights of legend and actual plot.
    When a numeric vector is supplied, the numbers are interpreted as
    \code{"lines"} units.  In addition, the unit or numeric vector
    may have named arguments (\samp{top} and \samp{bottom}), in which
    case the non-named argument specify the default value (recycled as
    needed), overloaded by the named arguments.}
  \item{legend_width}{An object of class \code{"unit"} of length
    1 specifying the width of the legend (if any). Default: 5 lines.}
  \item{main_gp, sub_gp}{object of class \code{"gpar"} containing the graphical parameters used for the main (sub) title, if specified.}

  \item{name}{name of the plotting viewport.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be fixed or not.}
  \item{prefix}{optional character string used as a prefix for the
    generated viewport and grob names.}
  \item{newpage}{logical. Should \code{\link{grid.newpage}} be called 
    before plotting?}
  \item{pop}{logical. Should the viewport created be popped?}
  \item{return_grob}{logical. Should a snapshot of the display be
    returned as a grid grob?}
  \item{\dots}{further arguments passed to \code{\link{grid_barplot}}.}
}
\details{
  The observed frequencies are displayed as bars and the fitted
  frequencies as a line.  By default a sqrt scale is used to make the
  smaller frequencies more visible.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org},
  David Meyer \email{David.Meyer@R-project.org}
}
\references{
  J. W. Tukey (1977),
  \emph{Exploratory Data Analysis}.
  Addison Wesley, Reading, MA.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\seealso{\code{\link{grid_barplot}}}
\examples{
## Simulated data examples:
dummy <- rnbinom(200, size = 1.5, prob = 0.8)
observed <- table(dummy)
fitted1 <- dnbinom(as.numeric(names(observed)),
                   size = 1.5, prob = 0.8) * sum(observed)
fitted2 <- dnbinom(as.numeric(names(observed)),
                   size = 2, prob = 0.6) * sum(observed)
rootogram(observed, fitted1)
rootogram(observed, fitted2)

## Real data examples:
data("HorseKicks")
HK.fit <- goodfit(HorseKicks)
summary(HK.fit)
plot(HK.fit)
## or equivalently
rootogram(HK.fit)

data("Federalist")
F.fit <- goodfit(Federalist, type = "nbinomial")
summary(F.fit)
plot(F.fit)

## (Pearson) residual-based shading
data("Federalist")
Fed_fit0 <- goodfit(Federalist, type = "poisson")
plot(Fed_fit0, shade = TRUE)
}
\keyword{hplot}
