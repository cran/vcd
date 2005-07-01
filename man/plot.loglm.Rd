\name{plot.loglm}
\alias{plot.loglm}
\alias{assoc.loglm}
\alias{mosaic.loglm}
\title{Visualize Fitted Log-linear Models}
\description{
  Visualize fitted \code{"loglm"} objects by mosaic or 
  association plots.
}
\usage{
\method{plot}{loglm}(x, panel = mosaic, type = c("observed", "expected"),
  residuals_type = c("pearson", "deviance"), gp = shading_hcl, gp_args = list(),
  \dots)
}
\arguments{
  \item{x}{a fitted \code{"loglm"} object, see \code{\link{loglm}}.}
  \item{panel}{a panel function for visualizing the observed values,
    residuals and expected values. Currently, \code{\link{mosaic}} and
    \code{\link{assoc}} in \pkg{vcd}.}
  \item{type}{a character string indicating whether the observed or the
    expected values of the table should be visualized.}
  \item{residuals_type}{a character string indicating the type of
    residuals to be computed.}
  \item{gp}{object of class \code{"gpar"}, shading function or a
    corresponding generating function (see details and
    \code{\link{shadings}}).
    Ignored if \code{shade = FALSE}.} 
  \item{gp_args}{list of arguments for the shading-generating function,
    if specified.}
  \item{\dots}{Other arguments passed to the \code{panel} function.}
}
\details{
  The \code{plot} method for \code{"loglm"} objects by default visualizes 
  the model using a mosaic plot (can be changed to an association plot
  by setting \code{panel = assoc}) with a shading based on the residuals
  of this model. The legend also reports the corresponding p value of the
  associated goodness-of-fit test. The \code{mosaic} and \code{assoc} methods
  are simple convenience interfaces to this \code{plot} method, setting 
  the \code{panel} argument accordingly.
}
\value{
  The \code{"structable"} visualized is returned invisibly.
}
\author{
  Achim Zeileis \email{Achim.Zeileis@R-project.org}
}
\seealso{
  \code{\link{loglm}},
  \code{\link{assoc}},
  \code{\link{mosaic}},
  \code{\link{strucplot}}
}
\examples{
## mosaic display for PreSex model
data(PreSex)
fm <- loglm(~ PremaritalSex * ExtramaritalSex * (Gender + MaritalStatus),
  data = aperm(PreSex, c(3, 2, 4, 1)))
fm
## visualize Pearson statistic
plot(fm, split_vertical = TRUE)
## visualize LR statistic
plot(fm, split_vertical = TRUE, residuals_type = "deviance")

## conditional independence in UCB admissions data
data(UCBAdmissions)
fm <- loglm(~ Dept * (Gender + Admit), data = aperm(UCBAdmissions))

## use mosaic display
plot(fm, labeling_args = list(abbreviate = c(Admit = 3)))

## and association plot
plot(fm, panel = assoc)
assoc(fm)
}
\keyword{hplot}
