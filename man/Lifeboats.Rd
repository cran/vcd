\name{Lifeboats}
\alias{Lifeboats}
\docType{data}
\title{Lifeboats on the Titanic}
\description{
Data from Mersey (1912) about the 18 (out of 20) lifeboats launched
before the sinking of the S. S. Titanic.
}

\usage{data(Lifeboats)}

\format{A data frame with 18 observations and 8 variables.
  \describe{
    \item{launch}{launch time in \code{"POSIXt"} format,}
    \item{side}{factor. Side of the boat,}
    \item{boat}{factor indicating the boat,}
    \item{crew}{number of male crew members on board,}
    \item{men}{number of men on board,}
    \item{women}{number of women (including female crew) on board,}
    \item{total}{total number of passengers,}
    \item{cap}{capacity of the boat.}
  }
}

\references{
L. Mersey (1912), Report on the loss of the ``Titanic'' (S. S.).
Parliamentary command paper 6452.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data:
\url{http://hotspur.psych.yorku.ca/ftp/sas/vcd/catdata/lifeboat.sas}}

\examples{
data(Lifeboats)
attach(Lifeboats)
ternaryplot(
  Lifeboats[,4:6],
  pch = ifelse(side=="Port", 1, 19),
  col = ifelse(side=="Port", "red", "blue"),
  id  = ifelse(men/total > 0.1, as.character(boat), NA),
  prop.size = 2,
  dimnames.position = "edge",
  main = "Lifeboats on the Titanic"
)
legend(
  0.7, 0.8,
  legend = c("SIDE", "Port", "Starboard"),
  pch = c(NA, 1, 19),
  col = c("black", "red", "blue"),
)
}

\keyword{datasets}
