\name{mosaicpairs}
\alias{mosaicpairs}
\alias{mosaicpairs.default}
\alias{mosaicpairs.formula}
\title{Mosaic Matrices}
\description{Plots matrices of mosaics.}
\usage{
\method{mosaicpairs}{default}(x, main = deparse(substitute(x)),
            xlab = NULL, ylab = NULL, labels, \dots,
            type = c("pairwise", "total", "conditional", "joint"),
            shade = TRUE, oma = NULL, cex.labels = NULL, label.pos = 0.5,
	    font.labels = 1, gap = 1)
\method{mosaicpairs}{formula}(formula, data = NULL, \dots, main, subset)
}
\arguments{
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute.  The table is
    best created by the \code{table()} command.}
  \item{main}{character string for the mosaic title.}
  \item{xlab,ylab}{x- and y-axis labels used for the plot; by default,
    the first and second element of \code{names(dimnames(X))} (i.e., the
    name of the first and second variable in \code{X}).}
  \item{labels}{the names of the variables.}
  \item{\dots}{further parameters (\emph{currently not used}).}
  \item{type}{a character string indicating the type of mosaic matrix,
    see details for more information.}
  \item{shade}{a logical indicating whether to produce extended mosaic
    plots, or a numeric vector of at most 5 distinct positive numbers
    giving the absolute values of the cut points for the residuals.  By
    default, \code{shade} is \code{TRUE}.  Using \code{shade = TRUE} cuts absolute values at 2 and
    4.}
  \item{oma}{A vector of the form c(bottom, left, top, right) giving the
    size of the outer margins in lines of text.}
  \item{cex.labels, font.labels}{graphics parameters for the text panel.}
  \item{label.pos}{y position of labels in the text panel.}
  \item{gap}{distance between subplots, in margin lines.}
  \item{formula}{a formula, such as \code{y ~ x}.}
  \item{data}{a data.frame (or list), or a contingency table from which
    the variables in \code{formula} should be taken.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used for plotting.}
}
\details{
  \code{mosaicpairs()} plots a matrix of pairwise mosaic plots, four
  independence types are distinguished:
  \code{"pairwise"}, \code{"total"}, \code{"conditional"} and \code{"joint"}.
  The pairwise mosaic matrix shows bivariate marginal relations, collapsed over
  all other variables. The total independence mosaic matrix shows mosaic plots
  for mutual independence, i.e. for marginal and conditional independence
  among all pairs of variables. The conditional independence mosaic matrix shows
  mosaic plots for marginal independence given all other variables. The
  joint independence mosaic matrix shows mosaic plots for joint independence of all pairs of
  variables from the others.
}

\references{
M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}
\seealso{
  \code{\link{mosaicplot}},
  \code{\link{pairs}}.
}
\examples{
data(UCBAdmissions)
data(PreSex)

mosaicpairs(PreSex, shade = TRUE)
mosaicpairs(UCBAdmissions, shade = TRUE)
mosaicpairs(UCBAdmissions, shade = TRUE, type="conditional")
}
\keyword{hplot}



