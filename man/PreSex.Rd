\name{PreSex}
\alias{PreSex}
\docType{data}
\title{Pre-marital Sex and Divorce}
\description{
  Data from Thornes \& Collard (1979), reported in Gilbert (1981), on
  pre- and extra-marital sex and divorce.
}
\usage{
data(PreSex)
}
\format{
  A 4-dimensional array resulting from cross-tabulating 1036
  observations on 4 variables.  The variables and their levels are as
  follows:

  \tabular{rll}{
    No \tab Name \tab Levels \cr
    1 \tab MaritalStatus \tab Divorced, Married \cr
    2 \tab ExtramaritalSex \tab Yes, No \cr
    3 \tab PremaritalSex \tab Yes, No \cr
    4 \tab Gender \tab Women, Men
  }
}

\references{
  G. N. Gilbert (1981),
  \emph{Modelling Society: An Introduction to Loglinear Analysis for
    Social Researchers}.
  Allen and Unwin, London.

  B. Thornes \& J. Collard (1979),
  \emph{Who Divorces?}.
  Routledge \& Kegan, London.

  M. Friendly (2000),
  \emph{Visualizing Categorical Data}.
  SAS Institute, Cary, NC.
}

\source{
  Michael Friendly (2000),
  Visualizing Categorical Data:
  \url{http://hotspur.psych.yorku.ca/ftp/sas/vcd/catdata/marital.sas}
}

\examples{
data(PreSex)

## Mosaic display for Gender and Premarital Sexual Experience
## (Gender Pre)
mosaic(margin.table(PreSex, c(3,4)), 
                main = "Gender and Premarital Sex")

## (Gender Pre)(Extra)
mosaic(margin.table(PreSex, c(2,3,4)), 
       expected = ~Gender * PremaritalSex + ExtramaritalSex ,
	   main = "PreMaritalSex*Gender +Sex")

## (Gender Pre Extra)(Marital)
mosaic(PreSex,
       expected = ~Gender*PremaritalSex*ExtramaritalSex + MaritalStatus,
       main = "PreMarital*ExtraMarital + MaritalStatus")

## (GPE)(PEM)
mosaic(PreSex, 
       expected = ~ Gender * PremaritalSex * ExtramaritalSex
                    + MaritalStatus * PremaritalSex * ExtramaritalSex,
       main = "G*P*E + P*E*M")
}

\keyword{datasets}
