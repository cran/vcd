\name{spacings}
\alias{spacings}
\alias{spacing_highlighting}
\alias{spacing_equal}
\alias{spacing_dimequal}
\alias{spacing_increase}
\alias{spacing_conditional}
\title{Spacing-generating Functions}
\description{
  These functions generate spacing functions to be used with
  \code{\link{strucplot}} to obtain customized spaces between the
  elements of a strucplot.
}
\usage{
spacing_equal(sp = unit(0.5, "lines"))
spacing_dimequal(sp)
spacing_increase(start = unit(0.3, "lines"), rate = 1.5)
spacing_conditional(sp = unit(0.5, "lines"), start = unit(2, "lines"), rate = 1.8)
spacing_highlighting(start = unit(0.2, "lines"), rate = 1.5)
}
\arguments{
  \item{start}{object of class \code{"unit"} indicating the start value
    for increasing spacings.}
  \item{rate}{increase rate for spacings.}
  \item{sp}{object of class \code{"unit"} specifying a fixed spacing.}
}
\details{
  These generating functions return a function used by
  \code{\link{strucplot}} to generate appropriate spaces between tiles of
  a strucplot, using the \code{dimnames} information of the visualized
  table.

  \code{spacing_equal} allows to specify one fixed space for \emph{all}
  dimensions.

  \code{spacing_dimequal} allows to specify a fixed space for
  \emph{each} dimension.

  \code{spacing_increase} creates increasing spaces for all dimensions,
  based on a starting value and an increase rate.

  \code{spacing_conditional} combines \code{spacing_equal} and
  \code{spacing_increase} to create fixed spaces for conditioned
  dimensions, and increasing spaces for conditioning dimensions.

  \code{spacing_highlighting} is essentially \code{spacing_conditional} but with
  the space of the last dimension set to 0. With a corresponding color
  scheme, this gives the impression of the last class being
  \sQuote{highlighted} in the penultimate class (as, e.g., in
  \code{\link{doubledecker}} plots).
}
\value{
  A spacing function with arguments:
  \item{d}{\code{"dim"} attribute of a contingency table.}
  \item{condvars}{index vector of conditioning dimensions (currently only used by
    \code{spacing_conditional}).}
  This function computes a list of objects of class \code{"unit"}.
  Each list element contains the spacing information for the
  corresponding dimension of the table.  The length of the
  \code{"unit"} objects is \eqn{k-1}, \eqn{k} number of levels of the
  corresponding factor.
}
\author{
  David Meyer \email{David.Meyer@R-project.org}
}
\seealso{
  \code{\link{strucplot}},
  \code{\link{doubledecker}}
}
\examples{
strucplot(Titanic, spacing = spacing_increase(start = 0.5, rate = 1.5))
strucplot(Titanic, spacing = spacing_equal(1))
strucplot(Titanic, spacing = spacing_dimequal(1:4 / 4))
strucplot(Titanic, spacing = spacing_highlighting,
                   gp = gpar(fill = c("light gray","dark gray")))
data(PreSex)
strucplot(aperm(PreSex, c(1,4,2,3)), spacing = spacing_conditional,
          condvars = 2)  
}
\keyword{hplot}
