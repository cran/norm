\name{prelim.norm}
\alias{prelim.norm}
\title{
Preliminary manipulations for a matrix of incomplete 
continuous data.
}
\description{
Sorts rows of x by missingness patterns, and centers/scales
columns of x. Calculates various bookkeeping quantities needed
for input to other functions, such as \code{em.norm} and \code{da.norm}.
}
\usage{
prelim.norm(x)
}
\arguments{
\item{x}{
data matrix containing missing values. The rows of x 
correspond to observational units, and the columns to variables. 
Missing values are denoted by \code{NA}.
}}
\value{
a list of thirteen components that summarize various features of x
after the data have been centered, scaled, and sorted by missingness
patterns. Components that might be of interest to the user include:

\item{nmis}{
a vector of length ncol(x) containing the number of missing
values for each variable in x. This vector has names that correspond
to the column names of x, if any.
}
\item{r}{
matrix of response indicators showing the missing data patterns 
in x.  Dimension is (S,p) where S is the number of distinct 
missingness patterns in the rows of x, and p is the number of 
columns in x. Observed values are indicated by 1 and missing 
values by 0. The row names give the number of observations in
each pattern, and the column names correspond to the column names of
x.
}}
\references{
See Section 5.3.1 of Schafer (1996).
}
\examples{
data(mdata)
s <- prelim.norm(mdata)  #do preliminary manipulations 
s$nmis[s$co] #look at nmis 
s$r #look at missing data patterns
}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
