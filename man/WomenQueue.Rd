\name{WomenQueue}
\alias{WomenQueue}
\docType{data}
\title{Women in Queues}
\description{
Data from Jinkinson & Slater (1981) and Hoaglin & Tukey (1985)
reporting the frequency distribution of females in 100 queues of
length 10 in a London Underground station.
}

\usage{data(WomenQueue)}

\format{A 1-way table giving the number of women
        in 100 queues of length 10. The variable and its levels are
  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab nWomen \tab 0, 1, ..., 10 \cr
  }
}


\references{
D. C. Hoaglin & J. W. Tukey (1985), Checking the shape of discrete distributions.
In D. C. Hoaglin, F. Mosteller, J. W. Tukey (eds.), \emph{Exploring Data Tables,
Trends and Shapes}, chapter 9. John Wiley & Sons, New York.

R. A. Jinkinson & M. Slater (1981), Critical discussion of a graphical method
for identifying discrete distributions, \emph{The Statistician}, \bold{30},
239-248.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data, pp. 19-20.}

\examples{
data(WomenQueue)
gf <- goodfit(WomenQueue, type = "binomial")
summary(gf)
plot(gf)
}

\keyword{datasets}
