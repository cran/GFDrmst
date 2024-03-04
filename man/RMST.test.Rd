\name{RMST.test}
\alias{RMST.test}
\title{
 Function to perform multiple RMST-based tests
}
\description{
 Linear hypotheses of the restricted mean survival times (RMSTs) of k different groups can be tested simultaneously.
 Therefore, the multivariate distribution of local Wald-type test statistics is approximated by (1) estimating the covariance between the test statistics (\code{method = "asymptotic"}), (2) a groupwise bootstrap procedure (\code{method = "groupwise"}) or (3) a permutation approach (\code{method = "permutation"}), respectively. Hence, adjusted p-values can be obtained.
}
\usage{
RMST.test(
  time = NULL,
  status = NULL,
  group = NULL,
  formula = NULL,
  event = NULL,
  data = NULL,
  hyp_mat,
  hyp_vec = NULL,
  tau,
  method = c("groupwise", "permutation", "asymptotic"),
  stepwise = FALSE,
  alpha = 0.05,
  Nres = 4999,
  seed = 1
)
}
\arguments{
  \item{time}{
  A vector containing the observed event times. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{status}{
  A vector of the same length as \code{time} containing the corresponding censoring status indicator with values 0 = censored and 1 = uncensored. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{group}{
  A vector of the same length as \code{time} containing the corresponding group labels. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{formula}{
  A model \code{formula} object. The left hand side contains the time variable and the
right hand side contains the factor variables of interest. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{event}{
  The name of censoring status indicator with values 0 = censored and 1 = uncensored. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{data}{
  A data.frame or list containing the variables in formula and the censoring status indicator. Default option is \code{NULL}. Either \code{time, status, group} or \code{formula, event, data} needs to be specified.
}
\item{hyp_mat}{
  A list containing all the hypothesis matrices H for the multiple tests or one of the options \code{"Tukey", "Dunnett"} and \code{"center"} for Tukey's or Dunnett's contrasts or the centering matrix, respectively, or a matrix if only one hypothesis is of interest. For the permutation test, all matrices need to be contrast matrices.
}
\item{hyp_vec}{
  A list containing all the hypothesis vectors c for the multiple tests or a vector if only one hypothesis is of interest. By default (\code{NULL}), all hypothesis vectors are set to zero vectors of suitable length.
}
\item{tau}{
  A numeric value specifying the end of the relevant time window for the analysis. Default option is \code{NULL}.
}
\item{method}{
 One of the methods \code{"groupwise", "permutation"} and \code{"asymptotic"} that should be used for calculating the critical values. Default option is \code{"groupwise"}.
}
\item{stepwise}{
 A logical vector indicating whether a stepwise extension of the test should be performed. If \code{TRUE}, no confidence intervals can be computed for the linear combinations but it may be that more tests can reject. Default option is \code{FALSE}.
}
\item{alpha}{
  A numeric value specifying the global level of significance. Default option is \code{0.05}.
}
\item{Nres}{
  The number of random variables to approximate the limiting distribution. This is only used if at least one hypothesis matrix is not a row vector. The default option is \code{4999}.
}
\item{seed}{
  A single value, interpreted as an integer, for providing reproducibility of the results or \code{NULL} if reproducibility is not wanted. Default option is \code{1}.
}
}
\details{
The restricted mean survival time (RMST) of group \eqn{i} is defined as
\deqn{ \mu_i := \int\limits_0^{\tau} S_i(t) \mathrm{d}t  }
for all \eqn{i\in\{1,...,k\}}, where \eqn{S_i} denotes the survival function of group \eqn{i}.
Let
\deqn{\boldsymbol{\mu} := (\mu_1,...,\mu_k)^{\prime}}
be the vector of the RMSTs and \deqn{\widehat{\boldsymbol{\mu}} := (\widehat{\mu}_1,...,\widehat{\mu}_k)^{\prime}} be the vector of their estimators.
Let \eqn{\mathbf{H}_{\ell} \in {\rm I\!R}^{r_{\ell}\times k}} with \eqn{\mathrm{rank}(\mathbf{H}_{\ell}) >0} and \eqn{\mathbf{c}_{\ell} \in {\rm I\!R}^{r_{\ell}}} for all \eqn{\ell\in\{1,...,L\}}.
We are considering the multiple testing problem with null and alternative hypotheses
\deqn{
\mathcal{H}_{0,\ell}: \mathbf{H}_{\ell} \boldsymbol\mu = \mathbf{c}_{\ell} \quad \mathrm{vs.} \quad \mathcal{H}_{1,\ell}: \mathbf{H}_{\ell} \boldsymbol\mu \neq \mathbf{c}_{\ell}, \qquad \mathrm{for }\: \ell\in\{1,...,L\}.
}

If \code{stepwise = TRUE}, the closed testing procedure is applied. In this case, no confidence intervals can be computed for the linear combinations but it may be that more tests can reject.

Note that the calculations for the groupwise bootstrap and the permutation test may take a while, especially with stepwise extension.
}
\value{
A list of class \code{GFDrmst} containing the following components:
\item{method}{A character containing the method which has been used.}
\item{test_stat}{A numeric vector containing the calculated Wald-type test statistics for the local hypotheses.}
\item{p.value}{A numeric vector containing the adjusted p-values for the local hypotheses.}
\item{res}{A list containing the results of the multiple tests including the hypothesis matrices, estimators of the linear combinations of RMSTs, potentially confidence intervals for the linear combinations (if all matrices are row vectors and stepwise = \code{FALSE}), Wald-type test statistics, critical values and the test decisions.}
\item{alpha}{A numeric value containing the global level of significance.}
}
\references{
Munko, M., Ditzhaus, M., Dobler, D., Genuneit, J. (2024). RMST-based multiple contrast tests in general factorial designs. Statistics in Medicine, 1-18. \doi{10.1002/sim.10017}
}
\examples{
data(colonCS, package = "condSURV")

# multiple asymptotic tests
out <- RMST.test(formula = "Stime ~ rx",
                 event = "event",
                 data = colonCS,
                 hyp_mat = "Tukey",
                 tau = 3000,
                 method = "asymptotic")
summary(out)
plot(out)
## or, equivalently,
out <- RMST.test(time = colonCS$Stime,
                 status = colonCS$event,
                 group = colonCS$rx,
                 hyp_mat = "Tukey",
                 tau = 3000,
                 method = "asymptotic")
summary(out)
plot(out)

\donttest{
# multiple groupwise bootstrap tests
# this may take a few seconds
out_gw <- RMST.test(formula = "Stime ~ rx",
                    event = "event",
                    data = colonCS,
                    hyp_mat = "Tukey",
                    tau = 3000,
                    method = "groupwise")
summary(out_gw)
plot(out_gw)

# multiple permutation tests
# this may take a few seconds
out_perm <- RMST.test(formula = "Stime ~ rx",
                      event = "event",
                      data = colonCS,
                      hyp_mat = "Tukey",
                      tau = 3000,
                      method = "permutation")
summary(out_perm)
plot(out_perm)
}
}
