\name{loddsratio}
\alias{loddsratio}
\alias{oddsratio}
\alias{loddsratio.default}
\alias{loddsratio.formula}

\alias{coef.loddsratio}
\alias{confint.loddsratio}
\alias{dim.loddsratio}
\alias{dimnames.loddsratio}
\alias{print.loddsratio}
\alias{vcov.loddsratio}
\alias{as.matrix.loddsratio}
\alias{as.array.loddsratio}
\alias{aperm.loddsratio}
\alias{t.loddsratio}
\alias{as.data.frame.loddsratio}

\title{
Calculate Generalized Log Odds Ratios for Frequency Tables
}
\description{
Computes (log) odds ratios and their asymptotic variance covariance
matrix for R x C (x strata) tables.
Odds ratios are calculated for two array dimensions, separately for each level of
all stratifying dimensions. See Friendly et al. (2011) for a sketch of a
general theory.
}
\usage{

loddsratio(x, \dots)
\method{loddsratio}{default}(x, strata = NULL, log = TRUE, 
  ref = NULL, correct = any(x == 0L), \dots)

\method{loddsratio}{formula}(formula, data = NULL, \dots,
subset = NULL, na.action = NULL)

oddsratio(x, stratum = NULL, log = TRUE)

\method{coef}{loddsratio}(object, log = object$log, \dots)
\method{vcov}{loddsratio}(object, log = object$log, \dots)
\method{print}{loddsratio}(x, log = x$log, \dots)
\method{confint}{loddsratio}(object, parm, level = 0.95, log = object$log, \dots)

%as.array(x, \dots)
\method{as.array}{loddsratio}(x, log=x$log, \dots)
\method{t}{loddsratio}(x)
\method{aperm}{loddsratio}(a, perm, \dots)
}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{an object. For the default method a k-way matrix/table/array
    of frequencies. The number of margins has to be at least 2.}
  \item{strata, stratum}{Numeric or character indicating the margins of a
    $k$-way table \code{x} (with $k$ greater than 2) that should be employed
    as strata. By default all dimensions except the first two are used.}
  \item{ref}{numeric or character. Reference categories for the (non-stratum)
    row and column dimensions that should be employed for computing the
    odds ratios. By default, odds ratios for profile contrasts (or sequential
    contrasts, i.e., successive differences of adjacent categories) are
    used. See details below.}
  \item{formula}{a formula specifying the variables used to create a 
    contingency table from \code{data}. A conditioning
    formula can be specified; the conditioning variables will then be
    used as strata variables.}
  \item{data}{either a data frame, or an object of class \code{"table"}
    or \code{"ftable"}.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s. Ignored if \code{data} is a contingency table.}
  \item{log}{logical. Should the results be displayed on a log scale
    or not? All internal computations are always on the log-scale but the
    results are transformed by default if \code{log = TRUE}.}
  \item{correct}{logical or numeric. Should a continuity correction
    be applied before computing odds ratios?
    If \code{TRUE}, 0.5 is  added to  all cells; if numeric (or an
    array conforming to the data) that value is added to all cells.  By default,
    this not  employed unless there are any zero cells in  the table, but this
    correction is often recommended to reduce bias when some frequencies are small
    (Fleiss, 1981).}
  \item{a, object}{an object of class \code{loddsratio} as computed by
    \code{loddsratio}.}
  \item{perm}{numeric or character vector specifying a permutation of strata.}
  \item{\dots}{arguments passed to methods.}
  
  \item{parm}{a specification of which parameters are to be given confidence intervals, 
    either a vector of numbers or a vector of names. If missing, all parameters are considered.}
  
  \item{level}{the confidence level required for the \code{confint} method.}

}
\details{
For an R x C table, (log) odds ratios are formed for the set of (R-1) x (C-1)
2 x 2 tables, corresponding to some set of contrasts among the row and column
variables.  The \code{ref} argument allows these to be specified in a general
way.

\code{ref = NULL} (default) corresponds to \dQuote{profile contrasts} 
(or sequential contrasts or successive differences) for ordered categories,
i.e., R1--R2, R2--R3, R3--R4, etc., and similarly for the column categories.
These are sometimes called \dQuote{local odds ratios}.

\code{ref = 1} gives contrasts with the first category; \code{ref = dim(x)}
gives contrasts with the last category; \code{ref = c(2, 4)} or \code{ref = list(2, 4)}
corresponds to the reference being the second category in rows and
the fourth in columns.

Combinations like \code{ref = list(NULL, 3)} are also possible, as are character
vectors, e.g., \code{ref = c("foo", "bar")} also works ("foo" pertaining again to the
    row reference and "bar" to column reference).

Note that all such parameterizations are equivalent, in that one can derive all
other possible odds ratios from any non-redundant set, but the interpretation
of these values depends on the parameterization.

Note also that these reference level parameterizations only have meaning when the
primary (non-strata) table dimensions are larger than 2x2.  In the 2x2 case,
the odds ratios are defined by the order of levels of those variables in the table,
so you can achieve a desired interpretation by manipulating the table.

See the help page of \code{\link{plot.loddsratio}} for visualization methods.
}
\value{
An object of class \code{loddsratio}, with the following components: 
\item{coefficients}{A named vector, of length (R-1) x (C-1) x \code{prod(dim(x)[strata])}
	containing the log odds ratios. Use the \code{coef} method to 
	extract these from the object, and the \code{confint} method for confidence intervals.
	For a two-way table, the names for the log odds ratios are constructed in the form
	Ri:Rj/Ci:Cj using the table names for rows and columns. For a stratified table,
	the names are constructed in the form Ri:Rj/Ci:Cj|Lk.
	}
\item{vcov}{Variance covariance matrix of the log odds ratios.}
\item{dimnames}{Dimension names for the log odds ratios, considered as a table of
  size (R-1, C-1, \code{dim(x)[strata]}). Use the \code{dim} and \code{dimnames} methods
  to extract these and manipulate the log odds ratios in relation to the original table.}
\item{dim}{Corresponding dimension vector.}
\item{contrasts}{A matrix C, such that \code{C \%*\% as.vector(log(x))} gives the log odds
  ratios. Each row corresponds to one log odds ratio, and is all zero, except for 4 elements
  of \code{c(1, -1, -1, 1)} for a given 2 x 2 subtable.}
\item{log}{A logical, indicating the value of \code{log} in the original call.}
}

\references{
A. Agresti (2013), \emph{Categorical Data Analysis}, 3rd Ed. New York: Wiley.

Fleiss, J. L. (1981). \emph{Statistical Methods for Rates and Proportions}.
2nd Edition. New York: Wiley.

M. Friendly (2000), \emph{Visualizing Categorical Data}. SAS Institute, Cary, NC.

Friendly, M., Turner, H,, Firth, D., Zeileis, A. (2011).
\emph{Advances in Visualizing Categorical Data Using the vcd, gnm and vcdExtra Packages in R}.
Correspondence Analysis and Related Methods (CARME 2011).
\url{http://www.datavis.ca/papers/adv-vcd-4up.pdf}

}
\author{
Achim Zeileis, Michael Friendly and David Meyer.
}
\note{
The method of calculation is an example of the use of the delta method
described by Agresti (2013), Section 16.1.6, giving estimates of
log odds ratios and their asymptotic covariance matrix.

The \code{coef} method returns the \code{coefficients} component as a vector
of length (R-1) x (C-1) x \code{prod(dim(x)[strata])}.
The \code{dim} and \code{dimnames} methods provide the proper attributes for
treating the \code{coefficients} vector as an (R-1) x (C-1) x strata array.
\code{as.matrix} and \code{as.array} methods are also provided for this purpose.

The \code{confint} method computes confidence intervals for the log odds ratios
(or for odds ratios, with \code{log = FALSE}).
The \code{\link[lmtest]{coeftest}} method (\code{summary} is an alias)
prints the asymptotic standard errors, z tests (standardized log odds
ratios), and the corresponding p values.

\emph{Structural zeros}: In addition to the options for zero cells provided by \code{correct},
the function allows for structural zeros to be represented as \code{NA} in the data argument.
\code{NA} in the data yields \code{NA} as the \code{LOR} estimate, but does not affect other
cells.

\code{oddsratio} is just an alias to \code{loddsratio} for backward
compatibility.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%\code{\link[vcd]{oddsratio}},
\code{\link{plot.loddsratio}} for some plotting methods;
\code{\link[stats]{confint}} for confidence intervals;
\code{\link[lmtest]{coeftest}} for z-tests of significance

}
\examples{
## artificial example
set.seed(1)
x <- matrix(rpois(5 * 3, 7), ncol = 5, nrow = 3)
dimnames(x) <- list(Row = head(letters, 3), Col = tail(letters, 5))

x_lor <- loddsratio(x)
coef(x_lor)
x_lor
confint(x_lor)
summary(x_lor)

## 2 x 2 x k cases
#data(CoalMiners, package = "vcd")
lor_CM <- loddsratio(CoalMiners)
lor_CM
coef(lor_CM)
confint(lor_CM)
confint(lor_CM, log = FALSE)

## 2 x k x 2
lor_Emp <-loddsratio(Employment)
lor_Emp
confint(lor_Emp)

## 4 way tables 
data(Punishment, package = "vcd")
lor_pun <- loddsratio(Freq ~ memory + attitude | age + education, data = Punishment)
lor_pun
confint(lor_pun)
summary(lor_pun)

# fit linear model using WLS
lor_pun_df <- as.data.frame(lor_pun)
pun_mod1 <- lm(LOR ~ as.numeric(age) * as.numeric(education),
               data = lor_pun_df, weights = 1 / ASE^2)
anova(pun_mod1)

## illustrate ref levels
VA.fem <- xtabs(Freq ~ left + right, subset=gender=="female", data=VisualAcuity)
VA.fem
loddsratio(VA.fem)                  # profile contrasts
loddsratio(VA.fem, ref=1)           # contrasts against level 1
loddsratio(VA.fem, ref=dim(VA.fem)) # contrasts against level 4


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{category}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
