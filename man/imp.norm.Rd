\name{imp.norm}
\alias{imp.norm}
\title{
Impute missing multivariate normal data
}
\description{
Draws missing elements of a data matrix under the multivariate normal
model and a user-supplied parameter
}
\usage{
imp.norm(s, theta, x)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix \code{x} created by the
function \code{prelim.norm}.
}
\item{theta}{
value of the normal parameter under which the missing data are to be
randomly imputed. This is a parameter vector in packed storage, such
as one created by \code{em.norm} or \code{da.norm}.
}
\item{x}{
the original data matrix used to create the summary list \code{s}. If this
argument is not supplied, then the data matrix returned by this
function may disagree slightly with the observed values in \code{x} due to
rounding errors.
}}
\value{
a matrix of the same form as \code{x}, but with all missing values filled in
with simulated values drawn from their predictive distribution given
the observed data and the specified parameter.
}
\details{
This function simply performs one I-step of data augmentation.
}
\section{WARNING}{
Before this function may be used, the random number generator seed
must be initialized with \code{rngseed} at least once in the current S
session.
}
\references{
See Section 5.4.1 of Schafer (1996).
}
\seealso{
\code{\link{prelim.norm}}, \code{\link{makeparam.norm}}, and \code{\link{rngseed}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
thetahat <- em.norm(s)   #find the mle
rngseed(1234567)   #set random number generator seed
ximp <- imp.norm(s,thetahat,mdata)  #impute missing data under the MLE
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
