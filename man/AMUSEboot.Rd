\name{AMUSEboot}
\alias{AMUSEboot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Second-order Separation Sub-White-Noise Bootstrap Testing with AMUSE
}
\description{
The function uses AMUSE to test, assuming a p-variate second-order stationary BSS model, whether the last \code{p-k} latent series are pure white noise. Four different bootstrapping strategies are available and the function can be run in parallel.
}
\usage{
AMUSEboot(X, k, tau = 1, n.boot = 200, s.boot = "p", ncores = NULL, iseed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A numeric data matrix or a numeric multivariate time series.}
  \item{k}{The number of latent series that are not white noise. Can be between 0 and p-1.}
  \item{tau}{The lag for the AMUSE autocovariance matrix.}
  \item{n.boot}{The number of bootstrapping samples.}
  \item{s.boot}{Bootstrapping strategy to be used. Possible values are \code{"p"}, \code{"np1"}, \code{"np2"}, \code{"np3"}. See details for further information.}
  \item{ncores}{The number of cores to be used. If \code{NULL} or 1, no parallel computing is used. Otherwise \link[parallel]{makeCluster} with \code{type = "PSOCK"} is used.}
  \item{iseed}{If parallel computation is used, the seed passed on to \code{\link[parallel]{clusterSetRNGStream}}. Default is \code{NULL} which means no fixed seed is used.}
}

\details{
AMUSE standardizes the data matrix \code{X} with \eqn{n} samples and computes the eigenedcomposition of the autocovariance matrix of the standardized data for a chosen lag \code{tau}, yielding a transformation \eqn{W} giving the latent variables as \eqn{S = XW}. Assume, without loss of generality, that the latent components are ordered in decreasing order with respect to the squares of the corresponding eigenvalues of the autocovariance matrix. Under the null the final \eqn{p - k} eigenvalues equal zero, \eqn{\lambda_{p-k} = \cdots = \lambda_{p}}{\lambda_(p-k) = ... = \lambda_p = 0}, and their mean square \eqn{m} can be used as a test statistic in bootstrap-based inference on the true number of latent white noise series.


The function offers four different bootstrapping strategies for generating samples for which the null approximately holds, and they are all based on the following general formula:

\enumerate{
\item Decompose the AMUSE-estimated latent series \eqn{S} into the postulated signal \eqn{S_1} and white noise \eqn{S_2}.
\item Take \eqn{n} bootstrap samples \eqn{S_2^*} of \eqn{S_2}, see the different strategies below.
\item Recombine \eqn{S^*=(S_1, S_2^*)} and back-transform \eqn{X^*= S^* W^-1}.
\item Compute the test statistic based on \eqn{X^*}. 
\item Repeat the previous steps \code{n.boot} times.
}

The four different bootstrapping strategies are:

\enumerate{
\item \code{s.boot = "p"}: 
The first strategy is parametric and simply generates all boostrap samples independently and identically from the standard normal distribution.
\item \code{s.boot = "np1"}: 
The second strategy is non-parametric and pools all observed \eqn{n(p - k)} white noise observations together and draws the bootstrap samples from amongst them.
\item \code{s.boot = "np2"}: 
The third strategy is non-parametric and proceeds otherwise as the second strategy but acts component-wise. That is, for each of the \eqn{p - k} white noise series it pools the observed \eqn{n} white noise observations together and draws the bootstrap samples of that particular latent series from amongst them.
\item \code{s.boot = "np3"}: 
The third strategy is non-parametric and instead of drawing the samples univariately as in the second and third strategies, it proceeds by resampling \eqn{n} vectors of size \eqn{p - k} from amongst all the observed \eqn{n} white noise vectors.
}

The function can be run in parallel by setting \code{ncores} to the desired number of cores (should be less than the number of cores available - 1). When running code in parallel the standard random seed of R is overridden and if a random seed needs to be set it should be passed via the argument \code{iseed}. The argument \code{iseed} has no effect in case \code{ncores} equals 1 (the default value).
}

\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{The value of the test statistic.}
  \item{p.value}{The p-value of the test.}
  \item{parameter}{The number of bootstrap samples.}
  \item{alternative}{Character string specifying the alternative hypothesis.}
  \item{k}{The number of latent series that are not white noise used in the testing problem.}
  \item{W}{The transformation matrix to the latent series.}
  \item{S}{Multivariate time series with the centered source components.}
  \item{D}{The underlying eigenvalues of the autocovariance matrix.}
  \item{MU}{The location of the data which was substracted before calculating AMUSE.}
  \item{tau}{The used lag.}
  \item{method}{Character string indicating which test was performed.}
  \item{data.name}{Character string giving the name of the data.}
  \item{s.boot}{Character string denoting which bootstrapping test version was used.}
}

\references{
Matilainen, M., Nordhausen, K., Virta, J. (2018), \emph{On the number of signals in multivariate time series}, to appear in the proceedings of International Conference on Latent Variable Analysis and Signal Separation 2018 (LVA/ICA), <https://arxiv.org/abs/1801.04925>.
}
\author{
Markus Matilainen, Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link[JADE]{AMUSE}}, \code{\link[JADE]{SOBI}}, \code{\link{SOBIboot}}
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

  boot_res_1 <- AMUSEboot(X, k = 1)
  boot_res_1

  boot_res_2 <- AMUSEboot(X, k = 2)
  boot_res_2

  # Plots of the estimated sources, the last two are white noise
  plot(boot_res_2$S)

}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ ts }
