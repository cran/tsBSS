\name{FixNA}
\encoding{latin1}
\alias{FixNA}
\alias{FixNA.default}
\alias{FixNA.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Variant of SOBI for Blind Source Separation
}
\description{
The FixNA (Shi et al., 2009) and FixNA2 (Matilainen et al., 2016) methods for blind source separation problem. It is used in a case of time series with stochastic volatility.
These methods are alternatives to vSOBI method.
}
\usage{

FixNA(X, ...)

\method{FixNA}{default}(X, k = 1:12, eps = 1e-06, maxiter = 1000, G = "pow", method = "FixNA", ...)
\method{FixNA}{ts}(X, \dots)

}
\arguments{
  \item{X}{Numeric matrix or multivariate time series object of class \code{\link{ts}}. Missing values are not allowed.}
  \item{k}{Vector of lags. Lag can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{Maximum number of iterations.}
  \item{G}{Function \eqn{G(x)}. Choices are \code{pow} (default) and \code{lcosh}.}
  \item{method}{Method to be used. Choices are \code{FixNA} (default) and \code{FixNA2}.}
  \item{...}{Further arguments to be passed to or from methods.}
}
\details{

Assume that \eqn{{\bf Y}}{Y} has \eqn{p} columns and it is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X} - \frac{1}{T}\sum_{t=1}^T {\bf X}_{ti})}{Y = S^(-1/2)*(X - (1/T)*sum_t(X_(ti)))}, where \eqn{{\bf S}}{S} is a sample covariance matrix of \eqn{{\bf X}}{X}. The algorithm for method FixNA finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing
\deqn{{\bf D}_1({\bf U}) = \sum_{k = 1}^K {\bf D}_{1k}({\bf U})= \sum_{k = 1}^K \sum_{i = 1}^p  \frac{1}{T - k}\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_t) G({\bf u}_i' {\bf Y}_{t + k})]}{D_1(U) = sum_k(D_k1(U)) = sum_i(sum_k((1/(T - k))*sum_t[G(u_i' Y_t)*G(u_i' Y_(t + k))])) }
and the algorithm for method FixNA2
\deqn{{\bf D}_2({\bf U}) = \sum_{k = 1}^K {\bf D}_{2k}({\bf U}) }{D_2(U) =  sum_k(D_k2(U))}
\deqn{= \sum_{k = 1}^K \sum_{i = 1}^p\left|\frac{1}{T - k}\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_t) G({\bf u}_i' {\bf Y}_{t + k})] - \left(\frac{1}{T - k}\right)^2\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_t)]\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_{t + k})]\right|.}{D_2(U) =  sum_k(D_k2(U)) = sum_k(sum_i(|(1/(T - k))*sum_t[G(u_i' Y_t)*G(u_i' Y_(t + k))] - (1/(T - k))^2*sum_t[G(u_i' Y_t)]*sum_t[G(u_i' Y_(t + k))]|)),}
\if{html}{where \eqn{i = 1, \ldots, p}, \eqn{k = 1, \ldots, K} and \eqn{t = 1, \ldots, T}.} For function \eqn{G(x)} the choices are \eqn{x^2} and log(cosh(\eqn{x})).

The algorithm works iteratively starting with \code{diag(p)} as an initial value for an orthogonal matrix \eqn{{\bf U} = ({\bf u}_1, {\bf u}_2, \ldots, {\bf u}_p)'}{U = (u_1, u_2, \ldots, u_p)'}.

Matrix \eqn{{\bf T}_{mik}}{T_mik} is a partial derivative of \eqn{{\bf D}_{mk}({\bf U})}{D_mk(U)}, for \eqn{m = 1, 2}, with respect to \eqn{{\bf u}_i}{u_i}.
Then \eqn{{\bf T}_{mk} = ({\bf T}_{m1k}, \ldots, {\bf T}_{mpk})'}{T_mk = (T_m1k, \ldots, T_mpk)'}, where \eqn{p} is the number of columns in \eqn{{\bf Y}}{Y}, and \eqn{{\bf T}_m = \sum_{k = 1}^K {\bf T}_{mk}}{T = sum(T_mk)}.
The update for the orthogonal matrix \eqn{{\bf U}_{new} = ({\bf T}_m{\bf T}_m')^{-1/2}{\bf T}_m}{U.new = (T_m T_m')^(-1/2)*T_m} is calculated at each iteration step. The algorithm stops when
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
\enc{Hyvärinen}{Hyvarinen}, A. (2001), \emph{Blind source separation by nonstationarity of variance: A cumulant-based approach}, IEEE Transactions on Neural Networks, 12(6): 1471-1474.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2016), \emph{ICA and Stochastic Volatility Models}, Proceedings of the XI International Conference on Computer Data Analysis and Modeling (CDAM'2016), 30--37.

Shi, Z., Jiang, Z. and Zhou, F. (2009), \emph{ Blind source separation with nonlinear autocorrelation
and non-gaussianity}, Journal of Computational and Applied Mathematics, 223(1): 908-915.

}
\author{
Markus Matilainen
}
\seealso{
\code{\link[tsBSS]{vSOBI}}
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
X <- ts(cbind(s1, s2, s3) \%*\% t(A), end = c(2015, 338), frequency = 365.25)


res <- FixNA(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}

\keyword{ multivariate }
\keyword{ ts }
