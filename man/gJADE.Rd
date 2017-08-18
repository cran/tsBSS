\name{gJADE}
\alias{gJADE}
\alias{gJADE.default}
\alias{gJADE.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized JADE
}
\description{
The gJADE method for blind source separation problem. It is designed for time series with stochastic volatility.
The method is a generalization of JADE, which is a method for blind source separation problem using only marginal information.
}
\usage{
gJADE(X, \ldots)

\method{gJADE}{default}(X, k = 0:12, eps = 1e-06, maxiter = 100, method = "frjd", \ldots)
\method{gJADE}{ts}(X, \ldots)

}
\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{A vector of lags. It can be any non-negative integer, or a vector consisting of them. Default is \code{0:12}. If \eqn{k = 0}, this method reduces to \code{\link[JADE]{JADE}}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{method}{The method to use for the joint diagonalization. The options are \code{\link[JADE]{rjd}} and \code{\link[JADE]{frjd}}. Default is \code{\link[JADE]{frjd}}.}
  \item{\ldots}{Other arguments passed on to chosen joint diagonalization method.}
}
\details{
Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}. The matrix \eqn{{\bf \widehat{C}}^{ij}_k({\bf Y})}{C^ij_k(Y)} is of the form

\deqn{{\bf \widehat{C}}^{ij}_k({\bf Y}) = {\bf \widehat{B}}^{ij}_k({\bf Y}) - {\bf S}_k({\bf Y}) ({\bf E}^{ij} + {\bf E}^{ji}) {\bf S}_k({\bf Y})' - \textrm{trace}({\bf E}^{ij}) {\bf I}_p,}{C^ij_k(Y) = B^ij_k(Y) - S_k(Y) (E^ij + E^ji) S_k(Y)' - trace(E^ij)*I,}

for \eqn{i, j = 1, \ldots, p}, where \eqn{{\bf S}_k({\bf Y})}{S_k(Y)} is the lagged sample covariance matrix of \eqn{{\bf Y}}{Y} for lag \eqn{k = 1, \ldots, K}, \eqn{{\bf E}^{ij}}{E^ij} is a matrix where element \eqn{(i,j)} equals to 1 and all other elements are 0, \eqn{{\bf I}_p}{I} is an identity matrix of order \eqn{p} and \eqn{{\bf \widehat{B}}^{ij}_k({\bf Y})}{B^ij_k(Y)} is as in \code{\link{gFOBI}}.

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing 
\deqn{\sum_{i = 1}^p \sum_{j = 1}^p \sum_{k = 0}^K ||diag({\bf U \widehat{C}}^{ij}_k({\bf Y}) {\bf U}')||^2.}{sum_i(sum_j (sum_k (||diag(U C^ij_k(Y) U')||^2))),}
\if{html}{where \eqn{k = 1, \ldots, K}.}
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.

}
\value{
A list with class 'bss' containing the following components:
  \item{W }{The estimated unmixing matrix.}
  \item{k }{The vector of the used lags.}
  \item{S }{The etimated sources as time series object standardized to have mean 0 and unit variances.}
}
\references{
Cardoso, J.-F., Souloumiac, A., (1993). \emph{Blind Beamforming for Non-Gaussian Signals}, in: IEE-Proceedings-F, volume 140, pp. 362--370.

Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New Independent Component Analysis Tools for Time Series}, Statistics & Probability Letters, 105, 80--87.

}
\author{
Klaus Nordhausen, Markus Matilainen
}
\seealso{
\code{\link[JADE]{frjd}}, \code{\link[JADE]{JADE}}, \code{\link{gFOBI}}
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

res <- gJADE(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}
\keyword{ multivariate }
\keyword{ ts }
