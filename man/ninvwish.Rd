\name{ninvwish}
\alias{ninvwish}
\title{
Random normal-inverted Wishart variate
}
\description{
Simulates a value from a normal-inverted Wishart distribution. This
function may be useful for obtaining starting values of the
parameters of a  multivariate normal distribution for multiple
chains of data augmentation.
}
\usage{
ninvwish(s, params)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix produced by the
function \code{prelim.norm}.
}
\item{params}{
list of parameters of a normal-inverted Wishart distribution. In
order, the elements of the list are: tau (a scalar), m (a scalar), mu0
(a vector of length ncol(x)), and lambdainv (a matrix of dimension
c(ncol(x),ncol(x))). When using this function to create starting
values for data augmentation, mu0 and lambdainv should be chosen in
relation to the data matrix after the columns have been centered and
scaled to have mean zero and variance one.
}}
\value{
a vector in packed storage representing the simulated normal-inverted
Wishart variate. This vector has the same form as parameter vectors
produced by functions such as \code{em.norm} and \code{da.norm}, and may be
used directly as a starting value for these functions. This vector can
also be put into a more understandable format by \code{getparam.norm}.
}
\section{WARNING}{
Before this function may be used, the random number generator seed
must be initialized with \code{rngseed} at least once in the current S
session.
}
\references{
See Section 5.4.2 of Schafer (1996).
}
\seealso{
\code{\link{rngseed}}, \code{\link{getparam.norm}}, \code{\link{em.norm}} and \code{\link{da.norm}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
params <- list(1,.5,rep(0,ncol(mdata)), .5*diag(rep(1,ncol(mdata)))) # gives widely dispersed values
rngseed(1234567)
start <- ninvwish(s,params)  # draw a variate
thetahat <- em.norm(s,start=start) # run EM from this starting value
}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
