\name{getparam.norm}
\alias{getparam.norm}
\title{
Extract normal parameters from packed storage
}
\description{
Takes a parameter vector, such as one produced by em.norm or da.norm,
and returns a list of parameters on the original scale.
}
\usage{
getparam.norm(s, theta, corr=FALSE)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix created by the
function \code{prelim.norm}.
}
\item{theta}{
vector of normal parameters expressed on transformed scale in packed
storage, such as one produced by the function \code{em.norm}.
}
\item{corr}{
if \code{TRUE}, computes means, standard deviations, and a correlation
matrix. If \code{FALSE}, computes means and a covariance matrix.
}}
\value{
if \code{corr=FALSE}, a list containing the components \code{mu} and \code{sigma}; if
\code{corr=TRUE}, a list containing the components \code{mu}, \code{sdv}, and \code{r}.  The
components are:

\item{mu}{
vector of means. Elements are in the same order and on the same scale
as the columns of the original data matrix, and with names
corresponding to the column names of the original data matrix.
}
\item{sigma}{
matrix of variances and covariances.
}
\item{sdv}{
vector of standard deviations.
}
\item{r}{
matrix of correlations.
}}
\seealso{
\code{\link{prelim.norm}} and \code{\link{makeparam.norm}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
thetahat <- em.norm(s)   #compute MLE
getparam.norm(s,thetahat,corr=TRUE)$r #look at estimated correlations
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
