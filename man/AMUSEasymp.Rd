\name{AMUSEasymp}
\alias{AMUSEasymp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second-order Separation Sub-White-Noise Asymptotic Testing with AMUSE
}
\description{
The function uses AMUSE to test, assuming a p-variate second-order stationary BSS model, whether the last \code{p-k} latent series are pure white noise. The test is asymptotic.
}
\usage{
AMUSEasymp(X, k, tau = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A numeric data matrix or a numeric multivariate time series.}
  \item{k}{The number of latent series that are not white noise. Can be between 0 and p-1.}
  \item{tau}{The lag for the AMUSE autocovariance matrix.}
}

\details{
AMUSE standardizes the data matrix \code{X} with \eqn{n} samples and computes the eigenedcomposition of the autocovariance matrix of the standardized data for a chosen lag \code{tau}, yielding a transformation \eqn{W} giving the latent variables as \eqn{S = XW}. Assume, without loss of generality, that the latent components are ordered in decreasing order with respect to the squares of the corresponding eigenvalues of the autocovariance matrix. Under the null the final \eqn{p - k} eigenvalues equal zero, \eqn{\lambda_{p-k} = \cdots = \lambda_{p}}{\lambda_(p-k) = ... = \lambda_p = 0}, and their mean square \eqn{m} can be used as a test statistic in inference on the true number of latent white noise series. 

This function conducts the hypothesis test using the asymptotic null distribution of \eqn{m}, a chi-squared distribution with \eqn{(p - k)(p - k + 1)/2} degrees of freedom.
}

\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{The value of the test statistic.}
  \item{p.value}{The p-value of the test.}
  \item{parameter}{The degrees of freedom of the asymptotic null distribution.}
  \item{method}{Character string indicating which test was performed.}
  \item{data.name}{Character string giving the name of the data.}
  \item{alternative}{Character string specifying the alternative hypothesis.}
  \item{k}{The number of latent series that are not white noise used in the testing problem.}
  \item{W}{The transformation matrix to the latent series.}
  \item{S}{Multivariate time series with the centered source components.}
  \item{D}{The underlying eigenvalues of the autocovariance matrix.}
  \item{MU}{The location of the data which was substracted before calculating AMUSE.}
  \item{tau}{The used lag.}
}

\references{
Nordhausen, K., Virta, J. (2018), \emph{Asymptotic test for the signal dimension in blind source separation}, under prepration.
}
\author{
Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link[JADE]{AMUSE}}, \code{\link[JADE]{SOBI}}, \code{\link{SOBIasymp}}
}
\examples{
  n <- 1000

  A <- matrix(rnorm(16), 4, 4)
  s1 <- arima.sim(list(ar = c(0.3, 0.6)), n)
  s2 <- arima.sim(list(ma = c(-0.3, 0.3)), n)
  s3 <- rnorm(n)
  s4 <- rnorm(n)

  S <- cbind(s1, s2, s3, s4)
  X <- S \%*\% t(A)

  asymp_res_1 <- AMUSEasymp(X, k = 1)
  asymp_res_1

  asymp_res_2 <- AMUSEasymp(X, k = 2)
  asymp_res_2

  # Plots of the estimated sources, the last two are white noise
  plot(asymp_res_2$S)

}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ ts }
