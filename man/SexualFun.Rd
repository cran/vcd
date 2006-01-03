\name{SexualFun}
\alias{SexualFun}
\docType{data}
\title{Sex is Fun}
\description{
  Data from Hout et al. (1987) given by Agresti (1990) summarizing the
  responses of married couples to the questionnaire item: Sex is fun for
  me and my partner: (a) never or occasionally, (b) fairly often, (c)
  very often, (d) almost always.
}
\usage{
data("SexualFun")
}
\format{
  A 2-dimensional array resulting from cross-tabulating the ratings of
  91 married couples.  The variables and their levels are as follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab Husband \tab Never Fun, Fairly Often, Very Often, Always Fun \cr
    2 \tab Wife \tab Never Fun, Fairly Often, Very Often, Always Fun
  }
}

\references{
  A. Agresti (1990),
  \emph{Categorical Data Analysis}.
  Wiley-Interscience, New York.

  M. Hout, O. D. Duncan, M. E. Sobel (1987),
  Association and heterogeneity: Structural models of similarities and
  differences,
  \emph{Sociological Methodology}, \bold{17}, 145-184.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  M. Friendly (2000),
  Visualizing Categorical Data, page 91.
}

\examples{
data("SexualFun")

## Kappa statistics
Kappa(SexualFun)

## Agreement Chart
agreementplot(t(SexualFun), weights = 1)
## Partial Agreement Chart and B-Statistics
agreementplot(t(SexualFun),
              xlab = "Husband's Rating",
              ylab = "Wife's Rating",
              main = "Husband's and Wife's Sexual Fun")
}

\keyword{datasets}
