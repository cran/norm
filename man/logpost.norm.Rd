\name{logpost.norm}
\alias{logpost.norm}
\title{
Observed-data log-posterior for normal data
}
\description{
Evaluates the log of the observed-data posterior density at a
user-supplied value of the parameter. Assumes a normal-inverted
Wishart prior. This function is useful for monitoring the progress of
EM and data augmentation.
}
\usage{
logpost.norm(s, theta, prior)
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
\item{prior}{
optional prior distribution. This is a list containing the
hyperparameters of a normal-inverted Wishart distribution. In order,
the elements of the list are: tau (a scalar), m (a scalar), mu0 (a
vector of length \code{ncol(x)}, where \code{x} is the original matrix of
incomplete data), and lambdainv (a matrix of dimension
\code{c(ncol(x),ncol(x))}). The elements of mu0 and lambdainv apply to the
data after transformation, i.e. after the columns have been centered
and scaled to have mean zero and variance one. If no prior is
supplied, the default is the usual noninformative prior for a
multivariate normal model: tau=0, m=-1, mu0=arbitrary, and lambdainv =
matrix of zeros.
}}
\value{
value of the observed-data log-posterior density
}
\references{
See Section 5.3.5 of Schafer (1996)
}
\seealso{
\code{\link{prelim.norm}} and \code{\link{loglik.norm}}
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
prior <- list(0,.5,rep(0,ncol(mdata)),
   .5*diag(rep(1,ncol(mdata)))) #ridge prior with .5 df
thetahat <- em.norm(s,prior=prior)   #compute posterior mode
logpost.norm(s,thetahat,prior)  #log-posterior at mode
}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
