\name{doubledecker}
\alias{doubledecker.default}
\alias{doubledecker.formula}
\alias{doubledecker}
\title{Doubledecker Plot}
\description{
  This function creates a doubledecker plot visualizing a classification
  rule.
}
\usage{
\method{doubledecker}{formula}(formula, data = NULL, \dots, main = NULL)
\method{doubledecker}{default}(x, depvar = length(dim(x)),
  margins = c(1,4, length(dim(x)) + 1, 1),
  col = rainbow_hcl(tail(dim(x), 1), end = 260),
  labeling = labeling_doubledecker,
  spacing = spacing_highlighting,
  main = NULL, keep_aspect_ratio = FALSE, \dots)
}
\arguments{
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}. The dependent variable is used
    last for splitting.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{x}{a contingency table in array form, with optional category
    labels specified in the \code{dimnames(x)} attribute.}
  \item{depvar}{dimension index or character string specifying the
    dependent variable.  That will be sorted last in the table.}
  \item{margins}{margins of the plot.  Note that by default, all factor
    names (except the last one) and their levels are visualized \emph{as
      a block} under the plot.}
  \item{col}{vector of colors used for the tiles of the last variable.}
  \item{labeling}{labeling function or corresponding generating
    generating function (see \code{\link{strucplot}} for details).}
  \item{spacing}{spacing object, spacing function or corresponding
    generating function (see \code{\link{strucplot}} for details).}
  \item{main}{either a logical, or a character string used for plotting
    the main title.  If \code{main} is \code{TRUE}, the
    name of the \code{data} object is used.}
  \item{keep_aspect_ratio}{logical indicating whether the aspect ratio should be
    maintained or not.}
  \item{\dots}{Further parameters passed to \code{mosaic}.}
}
\details{
  Doubledecker plots visualize the the dependence of one
  categorical (typically binary) variable on further categorical
  variables.  Formally, they are mosaic plots with vertical splits for
  all dimensions (antecedents) except the last one, which represents the
  dependent variable (consequent).  The last variable is visualized by
  horizontal splits, no space between the tiles, and separate colors for
  the levels.
}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\references{
  H. Hoffmann (2001),
  Generalized odds ratios for visual modeling.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{10}, 4, 628--640.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}

\seealso{
  \code{\link{strucplot}},
  \code{\link{mosaic}}
}
\examples{
doubledecker(Titanic)
doubledecker(Titanic, depvar = "Survived")
doubledecker(Survived ~ ., data = Titanic)
}
\keyword{hplot}

