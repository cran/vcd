\name{spineplot}
\alias{spineplot}
\alias{spineplot.default}
\alias{spineplot.formula}
\title{Spineplots}
\description{
  Spine plots are special cases of mosaic plots, \code{spineplot} is a
  convenience interface to \code{\link{mosaic}}.
}
\usage{
spineplot(x, ...)
\method{spineplot}{default}(x, split_vertical = TRUE, ...)
\method{spineplot}{formula}(formula, data = list(), cut = NULL,  ...)
}
\arguments{
  \item{x}{an object, the default method covers all cases covered by
    \code{\link{mosaic}}.}
  \item{formula}{a \code{"formula"} of type \code{y ~ x} with a single dependent \code{"factor"}
    and a single explanatory variable. Alternatively, this can be a \code{"table"}
    object.}    
  \item{split_vertical}{logical. Contrary to the settings in \code{\link{mosaic}}, the
    default is to start with a vertical split.}
  \item{data}{an optional data frame (only used if \code{formula} is a formula).}
  \item{cut}{if \code{x} is a numeric variable, it is \code{\link{cut}} into a
    factor at these points, by default set to \code{fivenum(x)}.}
  \item{...}{additional arguments passed to \code{\link{mosaic}}.}
}
\details{
  \code{spineplot} is a simple convenience interface to \code{mosaic}
  which uses a different default for the splitting direction. It is aimed
  at simple explanatory models like \code{y ~ x} where \code{x} can be 
  either a factor or a numeric which is \code{cut} to a factor.
}
\value{
  The structable visualized is returned invisibly.
}
\seealso{
  \code{\link{mosaic}}
}
\examples{
## Arthritis data (dependence on a categorical variable)
data(Arthritis)
(spineplot(Improved ~ Treatment, data = Arthritis))

## Space shuttle data (dependence on a numerical variable)
data(SpaceShuttle)
(spineplot(Fail ~ Temperature, data = SpaceShuttle))
}
\keyword{htest}
