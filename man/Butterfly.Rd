\name{Butterfly}
\alias{Butterfly}
\docType{data}
\title{Butterfly Species in Malaya}
\description{
Data from Fisher et al. (1943) giving the number of tokens found
for each of 501 species of butterflies collected in Malaya.
}

\usage{data(Butterfly)}

\format{A 1-way table giving the number of tokens for
        501 species of butterflies. The variable and its levels are
  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab nTokens \tab 0, 1, ..., 24 \cr
  }
}

\references{
R. A. Fisher, A. S. Corbet, C. B. Williams (1943), The relation between
the number of species and the number of individuals, \emph{Journal
of Animal Ecology}, \bold{12}, 42-58.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data, pp. 21-22.}

\examples{
data(Butterfly)
Ord.plot(Butterfly)
}

\keyword{datasets}
