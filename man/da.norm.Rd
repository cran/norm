\name{da.norm}
\alias{da.norm}
\title{
Data augmentation for incomplete multivariate normal data
}
\description{
Data augmentation under a normal-inverted Wishart prior. If no prior
is specified by the user, the usual "noninformative" prior for the
multivariate normal distribution is used. This function simulates one
or more iterations of a single Markov chain. Each iteration consists
of a random imputation of the missing data given the observed data and
the current parameter value (I-step), followed by a draw from the
posterior distribution of the parameter given the observed data and
the imputed data (P-step).
}
\usage{
da.norm(s, start, prior, steps=1, showits=FALSE, return.ymis=FALSE)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix produced by the
function \code{prelim.norm}.
}
\item{start}{
starting value of the parameter.  This is a parameter vector in packed
storage, such as one created by the function \code{makeparam.norm}. One
obvious choice for a starting value is an ML estimate or posterior
mode produced by \code{em.norm.}
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
}
\item{steps}{
number of data augmentation iterations to be simulated.
}
\item{showits}{
if \code{TRUE}, reports the iterations so the user can monitor the progress
of the algorithm.
}
\item{return.ymis}{
if \code{TRUE}, returns the output of the last I-step (imputed values of
missing data) in addition to the output of the last P-step. These
imputed values are useful for forming Rao-Blackwellized estimates of
posterior summaries.
}}
\value{
if \code{return.ymis=FALSE}, returns a parameter vector, the result of the last
P-step. If the value of \code{steps} was large enough to guarantee
approximate stationarity, then this parameter can be regarded as a
proper draw from the observed-data posterior, independent of \code{start}.
If \code{return.ymis=TRUE}, then this function returns a list of the following
two components:

\item{parameter}{
a parameter vector, the result of the last P-step
}
\item{ymis}{
a vector of missing values, the result of the last I-step.  The length
of this vector is \code{sum(is.na(x))}, where x is the original data
matrix. The storage order is the same as that of \code{x[is.na(x)]}.
}}
\section{WARNING}{
Before this function may be used, the random number generator seed
must be initialized with \code{rngseed} at least once in the current S
session.
}
\references{
See Chapter 5 of Schafer (1996).
}
\seealso{
\code{\link{rngseed}}, \code{\link{em.norm}}, \code{\link{prelim.norm}}, and \code{\link{getparam.norm}}.
}
\examples{
data(mdata)
s  <-  prelim.norm(mdata)
thetahat <- em.norm(s)   #find the MLE for a starting value
rngseed(1234567)   #set random number generator seed
theta <- da.norm(s,thetahat,steps=20,showits=TRUE)  # take 20 steps
getparam.norm(s,theta) # look at result
}
\keyword{distribution}
% Converted by Sd2Rd version 0.3-3.
