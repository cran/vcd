\name{agreementplot}
\alias{agreementplot}
\alias{agreementplot.default}
\alias{agreementplot.formula}
\title{Bangdiwala's Observer Agreement Chart}
\description{
  Representation of a \eqn{k \times k}{k by k} confusion matrix,
  where the observed and expected diagonal elements are represented by
  superposed black and white rectangles, respectively.  The function
  also computes a statistic measuring the strength of agreement
  (relation of respective area sums).
}
\usage{
\method{agreementplot}{default}(x, reverse_y = TRUE, main = NULL,
              weights = c(1, 1 - 1/(ncol(x) - 1)^2), margins = par("mar"),
              newpage = TRUE, pop = TRUE, xlab = names(dimnames(x))[2],
              ylab = names(dimnames(x))[1],
              xlab_rot = 0, xlab_just = "center",
              ylab_rot = 90, ylab_just = "center",
              fill_col = function(j) gray((1 - (weights[j]) ^ 2) ^ 0.5),
              line_col = "red", xscale = TRUE, yscale = TRUE,
               ...)
\method{agreementplot}{formula}(formula, data = NULL, ..., subset)
}
\arguments{
  \item{x}{a confusion matrix, i.e., a table with equal-sized dimensions.}
  \item{reverse_y}{if \code{TRUE}, the y axis is reversed (i.e., the
    rectangles' positions correspond to the contingency table).}
  \item{main}{user-specified main title.}
  \item{weights}{vector of weights for successive larger observed areas,
    used in the agreement strength statistic, and also for the
    shading.  The first element should be 1.}
  \item{margins}{vector of margins (see \code{\link[graphics]{par}}).}
  \item{newpage}{logical; if \code{TRUE}, the plot is drawn on a new page.}
  \item{pop}{logical; if \code{TRUE}, all newly generated viewports are popped after plotting.}
  \item{xlab, ylab}{labels of x- and y-axis.}
  \item{xlab_rot, ylab_rot}{rotation angle for the category labels.}
  \item{xlab_just, ylab_just}{justification for the category labels.}
  \item{fill_col}{a function, giving the fill colors used for exact and partial agreement}
  \item{line_col}{color used for the diagonal reference line}
  \item{formula}{a formula, such as \code{y ~ x}.
    For details, see \code{\link{xtabs}}.}
  \item{data}{a data frame (or list), or a contingency table from which
    the variables in \code{formula} should be taken.}
  \item{subset}{an optional vector specifying a subset of the rows in
    the data frame to be used for plotting.}
  \item{xscale, yscale}{logicals indicating whether the marginals should
  be added on the x-axis/y-axis, respectively.}
  \item{\dots}{further graphics parameters (see \code{\link{par}}).}
}
\details{
  Weights can be specified to allow for partial agreement, taking into
  account contributions from off-diagonal cells.  A weight vector of
  length 1 means strict agreement only, each additional element
  increases the maximum number of disagreement steps.
}
\value{
  Invisibly returned, a list with components
  \item{Bangdiwala}{the unweighted agreement strength statistic.}
  \item{Bangdiwala_Weighted}{the weighted statistic.}
  \item{weights}{the weight vector used.}
}
\references{
  Michael Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\examples{
data("SexualFun")
agreementplot(t(SexualFun))

data("MSPatients")
## best visualized using a resized device, e.g. using:
## get(getOption("device"))(width = 12)
pushViewport(viewport(layout = grid.layout(ncol = 2)))
pushViewport(viewport(layout.pos.col = 1))
agreementplot(t(MSPatients[,,1]), main = "Winnipeg Patients",
              newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
agreementplot(t(MSPatients[,,2]), main = "New Orleans Patients",
              newpage = FALSE)
popViewport(2)
dev.off()
}
\keyword{category}
\keyword{hplot}
