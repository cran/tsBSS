\name{vSOBI}
\alias{vSOBI}
\alias{vSOBI.default}
\alias{vSOBI.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Variant of SOBI for Blind Source Separation
}
\description{
The vSOBI method for blind source separation problem. It is used in case of time series with stochastic volatility.
The method is a variant of SOBI, which is a method designed to separate ARMA sources.
}
\usage{

vSOBI(X, ...)

\method{vSOBI}{default}(X, k = 1:12, eps = 1e-06, maxiter = 1000, ...)
\method{vSOBI}{ts}(X, \dots)

}
\arguments{
  \item{X}{Numeric matrix or multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{Vector of lags. Lag can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{Maximum number of iterations.}
  \item{...}{Further arguments to be passed to or from methods.}
}
\details{

Assume that \eqn{{\bf Y}}{Y} has \eqn{p} columns and it is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X} - \frac{1}{T}\sum_{t=1}^T {\bf X}_{ti})}{Y = S^(-1/2)*(X - (1/T)*sum_t(X_(ti)))}, where \eqn{{\bf S}}{S} is a sample covariance matrix of \eqn{{\bf X}}{X}. The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing
\deqn{{\bf D}({\bf U}) = \sum_{i = 1}^p \sum_{k = 1}^K \left(\frac{1}{T - k}\sum_{t=1}^T[({\bf u}_i' {\bf Y}_t)^2 ({\bf u}_i' {\bf Y}_{t + k})^2 - 1]\right)^2.}{D(U) = sum_i(sum_k(((1/(T - k))*sum_t[(u_i' Y_t)^2*(u_i' Y_(t + k)^2 - 1])^2)),}
\if{html}{where \eqn{i = 1, \ldots, p}, \eqn{k = 1, \ldots, K} and \eqn{t = 1, \ldots, T}.}

The algorithm works iteratively starting with \code{diag(p)} as an initial value for an orthogonal matrix \eqn{{\bf U} = ({\bf u}_1, {\bf u}_2, \ldots, {\bf u}_p)'}{U = (u_1, u_2, \ldots, u_p)'}. Matrix \eqn{{\bf T}_{ik}}{T_ik} is a partial derivative of \eqn{{\bf D}_k({\bf U})}{D_k(U)} with respect to \eqn{{\bf u}_i}{u_i}, where
\deqn{{\bf D}_k({\bf U}) = \sum_{i = 1}^p \left(\frac{1}{T - k}\sum_{t=1}^T[({\bf u}_i' {\bf Y}_t)^2 ({\bf u}_i' {\bf Y}_{t + k})^2 - 1]\right)^2.}{D_k(U) = sum_i(((1/(T - k))*sum_t[(u_i' Y_t)^2*(u_i' Y_(t + k)^2 - 1])^2).}

Then \eqn{{\bf T}_k = ({\bf T}_{1k}, \ldots, {\bf T}_{pk})'}{T_k = (T_1k, \ldots, T_pk)'}, where \eqn{p} is the number of columns in \eqn{{\bf Y}}{Y}, and \eqn{{\bf T} = \sum_{k = 1}^K {\bf T}_k}{T = sum(T_k)}.
The update for the orthogonal matrix \eqn{{\bf U}_{new} = ({\bf TT}')^{-1/2}{\bf T}}{U.new = (TT')^(-1/2)*T} is calculated at each iteration step. The algorithm stops when
\deqn{||{\bf U}_{new} - {\bf U}_{old}||}{||U.new - U.old||}
is less than \code{eps}.
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.


}
\value{
A list with class 'bss' containing the following components:
  \item{W }{The estimated unmixing matrix.}
  \item{k }{The vector of the used lags.}
  \item{S }{Estimated sources as time series object standardized to have mean 0 and unit variances.}
}
\references{
Belouchrani, A., Abed-Meriam, K., Cardoso, J.F. and Moulines, R. (1997), \emph{A blind source separation technique using second-order statistics}, IEEE Transactions on Signal Processing, 434--444.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2016), \emph{A variant of SOBI to accommodate stochastic volatility}, Submitted.

}
\author{
Markus Matilainen
}
\seealso{
\code{\link[JADE]{SOBI}}
}
\examples{
library(stochvol)
n <- 10000
A <- matrix(rnorm(9), 3, 3)

# simulate SV models
s1 <- svsim(n, mu = -10, phi = 0.8, sigma = 0.1)$y
s2 <- svsim(n, mu = -10, phi = 0.9, sigma = 0.2)$y
s3 <- svsim(n, mu = -10, phi = 0.95, sigma = 0.4)$y

# create a daily time series
X <- ts(cbind(s1, s2, s3) \%*\% t(A), end=c(2015, 338), frequency=365.25)


res <- vSOBI(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}

\keyword{ multivariate }
\keyword{ ts }
