\name{em.norm}
\alias{em.norm}
\title{
EM algorithm for incomplete normal data
}
\description{
Performs maximum-likelihood estimation on the matrix of incomplete
data using the EM algorithm. Can also be used to find a posterior mode
under a normal-inverted Wishart prior supplied by the user.
}
\usage{
em.norm(s, start, showits=TRUE, maxits=1000, criterion=0.0001, prior)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix produced by the
function \code{prelim.norm}.
}
\item{start}{
optional starting value of the parameter. This is a parameter vector
in packed storage, such as one created by the function
\code{makeparam.norm}. If no starting value is supplied, \code{em.norm} chooses
its own starting value.
}
\item{showits}{
if \code{TRUE}, reports the iterations of EM so the user can monitor the
progress of the algorithm.
}
\item{maxits}{
maximum number of iterations performed. The algorithm will stop if the
parameter still has not converged after this many iterations.
}
\item{criterion}{
convergence criterion. The algorithm stops when
the maximum relative difference in all of the estimated means,
variances, or covariances from one iteration to the next
is less than or equal to this value.
}
\item{prior}{
optional prior distribution. This is a list containing the
hyperparameters of a normal-inverted Wishart distribution. In order,
the elements of the list are: tau (a scalar), m (a scalar), mu0 (a
vector of length ncol(x)), and lambdainv (a matrix of dimension
c(ncol(x),ncol(x))). The elements of mu0 ans lambdainv apply to 
the data after transformation, i.e. after the columns have been
centered and scaled to have mean zero and variance one. If no prior is
supplied, the default is a uniform prior, which results in 
maximum-likelihood estimation.
}}
\value{
a vector representing the maximum-likelihood estimates of the normal
parameters. This vector contains means, variances, and covariances on
the transformed scale in packed storage.  The parameter can be
transformed back to the original scale and put into a more
understandable format by the function \code{getparam.norm}.
}
\details{
The default starting value takes all means on the transformed scale to
be equal to zero, and covariance matrix on the transformed scale equal
to the identity. All important computations are carried out in double
precision, using the sweep operator.
}
\references{
See Section 5.3 of Schafer (1994).
}
\seealso{
\code{\link{prelim.norm}}, \code{\link{getparam.norm}}, and \code{\link{makeparam.norm}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
thetahat <- em.norm(s)   #compute mle
getparam.norm(s,thetahat,corr=TRUE)$r #look at estimated correlations
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
