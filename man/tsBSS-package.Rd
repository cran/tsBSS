\name{tsBSS-package}
\alias{tsBSS-package}
\docType{package}
\title{
Tools for Blind Source Separation and Supervised Dimension Reduction for Time Series
}
\description{
Different estimates are provided to solve
the blind source separation problem for multivariate time series with stochastic volatility
and supervised dimension reduction problem for multivariate time series.
}

\details{
\tabular{ll}{
Package: \tab tsBSS\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2017-06-28\cr
License: \tab GPL (>= 2)
}
This package contains functions for blind source separation (BSS) problem for multivariate time series with stochastic volatility, like for example GARCH models.
The main functions of the package for BSS problem are
\itemize{
\item{\code{\link{FixNA}}}  {Function to solve the BSS problem. Algorithm is an alternative to \code{\link{vSOBI}} algorithm to acommodate stochastic volatility.}
\item{\code{\link{gFOBI}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{FOBI}} designed for time series with stochastic volatility.}
\item{\code{\link{gJADE}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{JADE}} designed for time series with stochastic volatility.}
\item{\code{\link{vSOBI}}}  {Function to solve the BSS problem. Algorithm is a variant of \code{\link[JADE]{SOBI}} algorithm and an alternative to \code{\link{FixNA}} to acommodate stochastic volatility.}
}

The main function of the package for supervised dimension reduction is
\itemize{
\item{\code{\link{tssdr}}}  {Function for supervised dimension reduction for multivariate time series. Includes methods TSIR, TSAVE and TSSH.}
}

Methods for ARMA models, like AMUSE and SOBI, and some non-stationary BSS methods for time series are implemented in the \code{\link[JADE]{JADE-package}}. See function \code{\link[dr]{dr}} for methods for supervised dimension reduction for iid observations.

}

\author{
Markus Matilainen,  Christophe Croux, Jari Miettinen, Klaus Nordhausen, Hannu Oja, Sara Taskinen

Maintainer: Markus Matilainen <markus.matilainen@utu.fi> 
}
\references{
Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New independent component analysis tools for time series}, Statistics & Probability Letters, 105, 80--87.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2017), \emph{On Independent Component Analysis with Stochastic Volatility Models}, Austrian Journal of Statistics, 46(3--4), 57--66.

Matilainen M., Croux C., Nordhausen K. & Oja H. (2017), \emph{Supervised Dimension Reduction for Multivariate Time Series}, Econometrics and Statistics, x(x), xx--xx. \doi{10.1016/j.ecosta.2017.04.002}

Matilainen M., Croux C., Nordhausen K. & Oja H. \emph{Sliced Average Variance Estimation for multivariate Time Series}. Manuscript.

Shi, Z., Jiang, Z. and Zhou, F. (2009), \emph{ Blind Source Separation with Nonlinear Autocorrelation and Non-Gaussianity}, Journal of Computational and Applied Mathematics, 223(1): 908--915.

}

\keyword{package}
\keyword{multivariate}
\keyword{ts}
