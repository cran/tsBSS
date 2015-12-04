\name{gFOBI}
\alias{gFOBI}
\alias{gFOBI.default}
\alias{gFOBI.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized FOBI
}
\description{
The gFOBI method for blind source separation problem. It is used in case of time series with stochastic volatility.
The method is a generalization of FOBI, which is a method designed for iid data.
}
\usage{
gFOBI(X, ...)

\method{gFOBI}{default}(X, k = 0:12, eps = 1e-06, maxiter = 100, method = "frjd", ...)
\method{gFOBI}{ts}(X, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Numeric matrix or multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{Vector of lags. Lag can be any non-negative integer, or a vector consisting of them. Default is \code{0:12}. If \eqn{k = 0}, this method reduces to \code{\link[JADE]{FOBI}}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{Maximum number of iterations.}
  \item{method}{Method to use for the joint diagonalization, options are \code{\link[JADE]{rjd}} and \code{\link[JADE]{frjd}}. Default is \code{\link[JADE]{frjd}}.}
  \item{\dots}{Other arguments passed on to chosen joint diagonalization method.}
}
\details{
Assume that \eqn{{\bf Y}}{Y} has \eqn{p} columns and it is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X} - \frac{1}{T}\sum_{t=1}^T {\bf X}_{ti})}{Y = S^(-1/2)*(X - (1/T)*sum_t(X_(ti)))}, where \eqn{{\bf S}}{S} is a sample covariance matrix of \eqn{{\bf X}}{X}. Algorithm first calculates
\deqn{{\bf \widehat{B}}^{ij}_k({\bf Y}) =  \frac{1}{T - k} \sum_{t = 1}^T [{\bf Y}_{t + k} {\bf Y}_t' {\bf E}^{ij} {\bf Y}_t {\bf Y}_{t + k}']}{B^ij_k(Y) = (1/(T - k))*sum[Y_(t + k) Y_t'  E^ij Y_t Y_(t + k)'],}
\if{html}{where \eqn{t = 1, \ldots, T},} and then 
\deqn{{\bf \widehat{B}}_k({\bf Y}) = \sum_{i = 1}^p {\bf \widehat{B}}^{ii}_k({\bf Y}).}{B_k(Y) = sum(B^ii_k(Y)),}
\if{html}{for \eqn{i = 1, \ldots, p}.}

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing 
\deqn{\sum_{k = 0}^K ||\textrm{diag}({\bf U \widehat{B}}_k({\bf Y}) {\bf U}')||^2.}{sum(||diag(U B_k(Y) U')||^2).}
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.
}

\value{
A list with class 'bss' containing the following components:
  \item{W }{The estimated unmixing matrix.}
  \item{k }{The vector of the used lags.}
  \item{S }{Estimated sources as time series object standardized to have mean 0 and unit variances.}
}
\references{

Cardoso, J.-F., (1989), \emph{ Source separation using higher order moments}, in: Proceedings of the IEEE International Conference on Acoustics, Speech and Signal Processing, 2109--2112.

Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New independent component analysis tools for time series}, Statistics & Probability Letters, 105, 80--87.

}
\author{
Markus Matilainen, Klaus Nordhausen
}

\seealso{
\code{\link[JADE]{FOBI}}, \code{\link[JADE]{frjd}}
}
\examples{
library(stochvol)
n <- 10000
A <- matrix(rnorm(9), 3, 3)

# simulate SV models
s1 <- svsim(n, mu = -10, phi = 0.8, sigma = 0.1)$y
s2 <- svsim(n, mu = -10, phi = 0.9, sigma = 0.2)$y
s3 <- svsim(n, mu = -10, phi = 0.95, sigma = 0.4)$y

X <- cbind(s1, s2, s3) \%*\% t(A)

res <- gFOBI(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero

}
\keyword{ multivariate }
\keyword{ ts }
