\name{structable}
\alias{structable.default}
\alias{structable.formula}
\alias{structable}
\alias{Extract.structable}
\alias{t.structable}
\alias{is.structable}
\alias{as.matrix.structable}
\alias{as.vector.structable}
\alias{dim.structable}
\alias{dimnames.structable}
\alias{as.table.structable}
\title{Structured Contingency Tables}
\description{
  This function produces a \sQuote{flat} representation of a
  high-dimensional contingency table constructed by recursive
  splits (similar to the construction of mosaic displays).
}
\usage{
\method{structable}{formula}(formula, data,
direction = NULL, split_vertical = NULL, \dots, subset, na.action)
\method{structable}{default}(\dots, direction = NULL, split_vertical = FALSE)
}
\arguments{
  \item{formula}{a formula object with possibly both left and right hand
    sides specifying the column and row variables of the flat table.}
  \item{data}{a data frame, list or environment containing the variables
    to be cross-tabulated, or an object inheriting from class \code{table}.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.
    Ignored if \code{data} is a contingency table.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.
    Ignored if \code{data} is a contingency table}
  \item{\dots}{\R objects which can be interpreted as factors (including
    character strings), or a list (or data frame) whose components can
    be so interpreted, or a contingency table object of class
    \code{"table"} or \code{"ftable"}.}
  \item{split_vertical}{logical vector indicating, for each dimension,
    whether it should be split vertically or not (default:
    \code{FALSE}). Values are recycled as needed.
    If the argument is of length 1, the value is alternated
    for all dimensions. Ignored if \code{direction} is provided.}
  \item{direction}{character vector alternatively specifying the
    splitting direction (\code{"h"} for horizontal and \code{"v"} for
    vertical splits).  Values are recycled as needed.  If the argument
    is of length 1, the value is alternated for all dimensions.}
}
\details{
  This function produces textual representations of mosaic displays, and
  thus \sQuote{flat} contingency tables. The formula interface is quite
  similar to the one of \code{\link{ftable}}, but also accepts the
  \code{\link{mosaic}}-like formula interface (empty left-hand
  side).  Note that even if the \code{\link{ftable}} interface is used,
  the \code{split_vertical} or \code{direction} argument is needed to
  specify the \emph{order} of the horizontal and vertical splits.
  If pretabulated data with a \code{Freq} column is used, than the
  left-hand side should be left empty---the \code{Freq} column will be
  handled correctly.

  \code{"structable"} objects can be subsetted using the \code{[[}
  operator, using either level indices or names (see examples).
}
\value{
  An object of class \code{"structable"},
  inheriting from class \code{"ftable"}, with the splitting
  information (\code{"split_vertical"}) as additional attribute.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{strucplot}},
  \code{\link{mosaic}},
  \code{\link[stats]{ftable}}
}

\examples{
structable(Titanic)
structable(Titanic, split_vertical = c(TRUE, TRUE, FALSE, FALSE))
structable(Titanic, direction = c("h","h","v","v"))
structable(Sex + Class ~ Survived + Age, data = Titanic)

## subsetting of structable objects
(hec <- structable(aperm(HairEyeColor)))

## subset _dimensions_ sequentially
hec[1]
hec[2]
hec["Male"]
hec[c("Male","Brown")]

## subset _rows_ and/or _columns_ sequentially
hec["Male","Hazel"]
hec[c("Male","Brown"),]
hec[c("Male","Brown"),"Hazel"]

## a few other operations
t(hec)
dim(hec)
dimnames(hec)
as.matrix(hec)

as.vector(hec) ## computed on the _multiway_ table
as.vector(unclass(hec))
}
\keyword{hplot}
