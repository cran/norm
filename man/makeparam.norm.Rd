\name{makeparam.norm}
\alias{makeparam.norm}
\title{
Convert normal parameters to packed storage
}
\description{
Does the opposite of \code{getparam.norm}.
Converts a list of user-specified parameters to a parameter vector
suitable for input to functions such as \code{da.norm} and \code{em.norm}.
}
\usage{
makeparam.norm(s, thetalist)
}
\arguments{
\item{s}{
summary list of an incomplete normal data matrix created
by the function \code{prelim.norm}.
}
\item{thetalist}{
list of normal parameters of the same form as one produced by
\code{getparam.norm}. If the list has two components, the first must be the
vector of means and the second must be the covariance matrix,
where means and covariances are expressed on the scale of the original
data. If the list has three components, the first must be the vector
of means, the second must be the vector of standard deviations, and
the third must be the correlation matrix.
}}
\value{
normal parameter in packed storage, suitable for use as a starting
value for \code{em.norm}, \code{mda.norm}, or \code{mdamet.norm}.
}
\seealso{
\code{\link{prelim.norm}} and \code{\link{getparam.norm}}.
}
\examples{
data(mdata)
s <- prelim.norm(mdata)   #do preliminary manipulations
thetahat <- em.norm(s)   #compute mle
thetahat <- getparam.norm(s,thetahat,corr=TRUE)   #extract parameters
thetahat$r   #look at mle correlations
thetahat$r[1,2] <- .5   #tweak a parameter
thetahat <- makeparam.norm(s,thetahat)   #convert to packed storage
thetahat <- em.norm(s,thetahat) #run EM again from new starting value
}
\keyword{regression}
% Converted by Sd2Rd version 0.3-3.
