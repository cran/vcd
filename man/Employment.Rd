\name{Employment}
\alias{Employment}
\docType{data}
\title{Employment Status}
\description{
Data from a 1974 Danish study given by Anderson (1991) on the
employees who had been laid off. The workers are classified by
their employment status on 1975-01-01, the cause of their layoff
and the length of employment before they were laid off.
}

\usage{data(Employment)}

\format{A 3-dimensional array resulting from cross-tabulating
  variables for 1314 employees. The variables and their levels are as
  follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab EmploymentStatus \tab NewJob, Unemployed \cr
    2 \tab EmploymentLength \tab <1Mo, 1-3Mo, 3-12Mo, 1-2Yr, 2-5Yr, >5Yr \cr
    3 \tab LayoffCause \tab Closure, Replaced
  }
}

\references{
E. B. Anderson (1991), \emph{Statistical Analysis of Categorical Data}. Springer-Verlag, Berlin.

M. Friendly (2000),
\emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.
}

\source{Michael Friendly (2000), Visualizing Categorical Data, p. 126-129.
}

\examples{
data(Employment)

## Employment Status ##
mosaicplot(Employment, shade = TRUE,
           margin = ~ LayoffCause*EmploymentLength + EmploymentStatus,
	   main = "Layoff*EmployLength + EmployStatus")


mosaicplot(Employment, shade = TRUE,
           margin = ~ LayoffCause*EmploymentLength + LayoffCause*EmploymentStatus,
	   main = "Layoff*EmployLength + Layoff*EmployStatus")

par(mfrow=c(1,2))

## Closure ##
mosaicplot(Employment[,,1], shade = TRUE, main = "Layoff: Closure")

## Replaced ##
mosaicplot(Employment[,,2], shade = TRUE, main = "Layoff : Replaced")

par(mfrow=c(1,1))
}

\keyword{datasets}
