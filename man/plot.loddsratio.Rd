\name{plot.loddsratio}
\alias{plot.loddsratio}
\alias{lines.loddsratio}
\title{Plotting (Log) Odds Ratios}
\description{
  Produces a (conditional) line plot of extended (log) odds ratios.
}
\usage{
\method{plot}{loddsratio}(x, baseline = TRUE, gp_baseline = gpar(lty = 2),
     lines = TRUE, lwd_lines = 3,
     confidence = TRUE, conf_level = 0.95, lwd_confidence = 2,
     whiskers = 0, transpose = FALSE,
     col = NULL, cex = 0.8, pch = NULL,
     bars = NULL, gp_bars = gpar(fill = "lightgray", alpha = 0.5),
     bar_width = unit(0.05, "npc"),
             
     legend = TRUE, legend_pos = "topright", legend_inset = c(0, 0),
     legend_vgap = unit(0.5, "lines"),
     gp_legend_frame = gpar(lwd = 1, col = "black"),
     gp_legend_title = gpar(fontface = "bold"),
     gp_legend = gpar(), legend_lwd = 1, legend_size = 1,
             
     xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
     main = NULL, gp_main = gpar(fontsize = 12, fontface = "bold"),
     newpage = TRUE, pop = FALSE, return_grob = FALSE,
     add = FALSE, prefix = "", \dots)
\method{lines}{loddsratio}(x, legend = FALSE, confidence = FALSE, cex = 0, \dots)
}
\arguments{
  \item{x}{an object of class \code{loddsratio}.}
  \item{baseline}{if \code{TRUE}, a dashed line is plotted at a
    value of 1 (in case of odds) or 0 (in case of log-odds).}
  \item{gp_baseline}{object of class \code{"gpar"} used for the baseline.}
  \item{lines}{if \code{TRUE}, the points are connected by lines (only
    sensible if the variable represented by the x-axis is ordinal).}
  \item{lwd_lines}{Width of the connecting lines (in \code{char} units).}
  \item{confidence}{logical; shall confindence intervals be plotted?}
  \item{conf_level}{confidence level used for confidence intervals.}
  \item{lwd_confidence}{Line width of the confidence interval bars (in \code{char} units).}
  \item{whiskers}{width of the confidence interval whiskers.}
  \item{transpose}{if \code{TRUE}, the plot is transposed.}

  \item{col}{character vector specifying the colors of the fitted
    lines, by default chosen with \code{\link[colorspace]{rainbow_hcl}}.}
  \item{cex}{size of the plot symbols (in lines).}
  \item{pch}{character or numeric vector of symbols used for plotting
    the (possibly conditioned) observed values, recycled as needed.}
  \item{bars}{logical; shall bars be plotted additionally to the points?
    Defaults to \code{TRUE} in case of only one conditioning variable.}
  \item{gp_bars}{object of class \code{"gpar"} used for the bars.}
  \item{bar_width}{Width of the bars, if drawn.}
  
  \item{legend}{logical; if \code{TRUE} (default), a legend is drawn.}
  \item{legend_pos}{numeric vector of length 2, specifying x and y
    coordinates of the legend, or a character string (e.g., \code{"topleft"},
    \code{"center"} etc.). Defaults to \code{"topleft"}
    if the fitted curve's slope is
    positive, and \code{"topright"} else.}
  \item{legend_inset}{numeric vector or length 2 specifying the inset
    from the legend's x and y coordinates in npc units.}
  \item{legend_vgap}{vertical space between the legend's line entries.}
  \item{gp_legend_frame}{object of class \code{"gpar"} used for the
    legend frame.}
  \item{gp_legend_title}{object of class \code{"gpar"} used for the
    legend title.}
  \item{gp_legend}{object of class \code{"gpar"} used for the
    legend defaults.}
  \item{legend_lwd}{line width used in the legend for the different
    groups.}
  \item{legend_size}{size used for the group symbols (in char units).}

  \item{xlab}{label for the x-axis. Defaults to \code{"Strata"} if
    \code{transpose} is \code{FALSE}.}
  \item{ylab}{label for the y-axis. Defaults to \code{"Strata"} if
    \code{transpose} is \code{TRUE}.}
  \item{xlim}{x-axis limits. Ignored if \code{transpose} is \code{FALSE}.}
  \item{ylim}{y-axis limits. Ignored if \code{transpose} is \code{TRUE}.}
  \item{main}{user-specified main title.}
  \item{gp_main}{object of class \code{"gpar"} used for the main title.}
  \item{newpage}{logical; if \code{TRUE}, the plot is drawn on a new page.}
  \item{pop}{logical; if \code{TRUE}, all newly generated viewports are
    popped after plotting.}
  \item{return_grob}{logical. Should a snapshot of the display be
    returned as a grid grob?}
  \item{add}{logical; should the plot added to an existing log odds
    ratio plot?}
  \item{prefix}{character string used as prefix for the viewport name.}
  \item{\dots}{other graphics parameters (see \code{\link{par}}).}
}
\value{
  if \code{return_grob} is \code{TRUE}, a grob object corresponding to
  the plot. \code{NULL} (invisibly) else.
}
\details{
  The function basically produces conditioned line plots of the (log)
  odds ratios structure provided in \code{x}.

  The \code{lines} method can be used to overlay different plots (for
  example, observed and expected values).

  \code{\link{cotabplot}} can be used for stratified analyses (see examples).

  }
\references{
  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\seealso{
  \code{\link{loddsratio}}
}
\examples{
## 2 x 2 x k cases
data(CoalMiners, package = "vcd")
lor_CM <- loddsratio(CoalMiners)
plot(lor_CM)
lor_CM_df <- as.data.frame(lor_CM)

# fit linear models using WLS
age <- seq(20, 60, by = 5)
lmod <- lm(LOR ~ age, weights = 1 / ASE^2, data = lor_CM_df)
grid.lines(seq_along(age), fitted(lmod), gp = gpar(col = "blue", lwd = 2), default.units = "native")
qmod <- lm(LOR ~ poly(age, 2), weights = 1 / ASE^2, data = lor_CM_df)
grid.lines(seq_along(age), fitted(qmod), gp = gpar(col = "red", lwd = 2), default.units = "native")

## 2 x k x 2
lor_Emp <-loddsratio(Employment)
plot(lor_Emp)

## 4 way tables 
data(Punishment, package = "vcd")
mosaic(attitude ~ age + education + memory, data = Punishment,
highlighting_direction="left", rep = c(attitude = FALSE))

# visualize the log odds ratios, by education
plot(loddsratio(~ attitude + memory | education, data = Punishment))

# visualize the log odds ratios, by age
plot(loddsratio(~ attitude + memory | age, data = Punishment))

# visualize the log odds ratios, by age and education
plot(loddsratio(~ attitude + memory | age + education, data = Punishment))

# same, transposed
plot(loddsratio(~ attitude + memory | age + education, data = Punishment), transpose = TRUE)

# alternative visualization methods
image(loddsratio(Freq ~ ., data = Punishment))
tile(loddsratio(Freq ~ ., data = Punishment))

## cotabplots for more complex tables
cotabplot(Titanic, cond = c("Age","Sex"), panel = cotab_loddsratio)

cotabplot(Freq ~ opinion + grade + year | gender, data = JointSports, panel = cotab_loddsratio)
cotabplot(Freq ~ opinion + grade | year + gender, data = JointSports, panel = cotab_loddsratio)
}
\keyword{category}
