\name{rngseed}
\alias{rngseed}
\title{
Initialize random number generator seed
}
\description{
Initializes the seed value for the internal random-number generator
used in missing-data programs
}
\usage{
rngseed(seed)
}
\arguments{
\item{seed}{
a positive number > = 1, preferably a large integer.
}}
\value{
\code{NULL}.
}
\section{NOTE}{
The random number generator seed must be set at least once
by this function before the simulation or imputation functions
in this package (\code{da.norm}, etc.) can be used.
}
\keyword{
distribution
}
% Converted by Sd2Rd version 0.3-3.
