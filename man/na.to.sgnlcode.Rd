\name{.na.to.snglcode}
\alias{.na.to.snglcode}
\title{
Changes NA's to single precision missing value code
}
\description{
Changes NA's to single precision missing value code
It's called internally by other functions in the package, e.g., 
\emph{'prelim.norm'}.
}
\usage{
.na.to.snglcode(x, mvcode)
}
\arguments{
\item{x}{data object.
}
\item{mvcode}{internal input of 'prelim.norm'.
}
}
\value{
Initial data object with missing values code precision changed to
sinlge.
}
\seealso{
\code{\link{prelim.norm}}
}
\keyword{
NA
}
