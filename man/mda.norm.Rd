\name{mda.norm}
\alias{mda.norm}
\title{
Monotone data augmentation for incomplete multivariate normal data
}
\description{
Monotone data augmentation under the usual noninformative prior, as
described in Chapter 6 of Schafer (1996).  This function simulates one
or more iterations of a single Markov chain. One iteration consists of
a random imputation of the missing data given the observed data and
the current parameter value (I-step), followed by a draw from the
posterior distribution of the parameter given the observed data and
the imputed data (P-step). The I-step imputes only enough data to
complete a monotone pattern, which typically makes this algorithm
converge more quickly than \code{da.norm}, particularly when the observed
data are nearly monotone. The order of the variables in the original
data matrix determines the monotone pattern to be completed. For fast
convergence, it helps to order the variables according to their rates
of missingness, with the most observed (least missing) variable on the
left and the least observed variable on the right.
}
\usage{
mda.norm(s, theta, steps=1, showits=FALSE)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix produced by the
function \code{prelim.norm}.
}
\item{theta}{
starting value of the parameter.  This is a parameter vector in packed
storage, such as one created by the function \code{makeparam.norm}. One
obvious choice for a starting value is an ML estimate or posterior
mode produced by \code{em.norm.}
}
\item{steps}{
number of monotone data augmentation iterations to be simulated.
}
\item{showits}{
if \code{TRUE}, reports the iterations so the user can monitor the progress
of the algorithm.
}}
\value{
Returns a parameter vector, the result of the last P-step. If the
value of \code{steps} was large enough to guarantee approximate
stationarity, then this parameter can be regarded as a proper draw
from the observed-data posterior, independent of \code{start}.
}
\section{WARNING}{
Before this function may be used, the random number generator seed
must be initialized with \code{rngseed} at least once in the current S
session.
}
\references{
Chapter 6 of Schafer (1996).
}
\seealso{
\code{\link{rngseed}}, \code{\link{em.norm}}, \code{\link{prelim.norm}}, and \code{\link{getparam.norm}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)
thetahat <- em.norm(s)   #find the MLE for a starting value
rngseed(1234567)   #set random number generator seed
theta <- mda.norm(s,thetahat,steps=20,showits=TRUE)  # take 20 steps
getparam.norm(s,theta) # look at result
}
\keyword{multivariate}
% Converted by Sd2Rd version 0.3-3.
