\name{summary.GFDrmst}
\alias{summary.GFDrmst}
\title{
Summarizing RMST-based multiple testing output
}
\description{
\code{summary} method for class \code{"GFDrmst"}.
}
\usage{
\method{summary}{GFDrmst}(object, digits = 8,...)
}
\arguments{
  \item{object}{
 An object of class \code{"GFDrmst"}, usually a result of a call to \code{\link{RMST.test}}.
}
  \item{digits}{
An integer indicating the number of decimal places to be used. Default option is \code{8}.
}
\item{...}{Further arguments are ignored.}
}
\details{
The function prints the information about the used method, significance level, hypothesis matrices, Wald-type test statistics, adjusted p-values and the overall results of the tests.
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
summary(out, digits = 3)
}
