\name{plot.GFDrmst}
\alias{plot.GFDrmst}
\title{
Plot simultaneous confidence intervals for a \code{GFDrmst} object
}
\description{
Simultaneous confidence intervals for the linear combinations of restricted mean survival times are plotted.
}
\usage{
\method{plot}{GFDrmst}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"GFDrmst"}, usually a result of a call to \code{\link{RMST.test}}.
}
\item{...}{Further arguments are ignored.}
}
\details{
The displayed vectors on the y-axis are the coefficients \eqn{\mathbf{H}_{\ell}} for the linear combinations of the RMSTs.
}
\value{
No return value, called for side effects.
}

\seealso{
\code{\link{RMST.test}}
}
\examples{
data(colonCS, package = "condSURV")

out <- RMST.test(formula = "Stime ~ rx",
                 event = "event",
                 data = colonCS,
                 hyp_mat = "Tukey",
                 tau = 3000,
                 method = "asymptotic")
plot(out)
}
