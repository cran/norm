\name{.code.to.na}
\alias{.code.to.na}
\title{
Changes missing value code to NA
}
\description{
Changes missing value code to NA. It's called from
\emph{'prelim.norm'}.
}
\usage{
.code.to.na(x, mvcode)
}
\arguments{
\item{x}{data object.
}
\item{mvcode}{internal input of 'prelim.norm'.
}
}
\value{
Initial data object with missing values code changed to NA.
}
\seealso{
\code{\link{prelim.norm}}
}
\keyword{
NA
}
