\name{loglik.norm}
\alias{loglik.norm}
\title{
Observed-data loglikelihood for normal data
}
\description{
Evaluates the observed-data loglikelihood function at a user-supplied
value of the parameter. This function is useful for monitoring the
progress of EM and data augmentation.
}
\usage{
loglik.norm(s, theta)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix created by the
function \code{prelim.norm}.
}
\item{theta}{
vector of normal parameters expressed on transformed scale in packed
storage, such as one produced by the function \code{em.norm}.
}}
\value{
value of the observed-data loglikelihood 
}
\references{
See Section 5.3.5 of Schafer (1996)
}
\seealso{
\code{\link{prelim.norm}} and \code{\link{logpost.norm}}
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
thetahat <- em.norm(s)   #compute MLE
loglik.norm(s,thetahat)  #loglikelihood at the MLE
}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
