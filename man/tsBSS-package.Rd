\name{tsBSS-package}
\alias{tsBSS-package}
\docType{package}
\title{
Tools for Blind Source Separation for Time Series
}
\description{
Different estimates are provided to solve the blind source separation problem for time series with stochastic volatility.
}

\details{
\tabular{ll}{
Package: \tab tsBSS\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2016-09-01\cr
License: \tab GPL (>= 2)
}
This package contains functions for blind source separation (BSS) problem for multivariate time series with stochastic volatility, like for example GARCH models.
The main functions of the package are
\itemize{
\item{\code{\link{FixNA}}}  {Function to solve the BSS problem. Algorithm is an alternative to \code{\link{vSOBI}} algorithm to acommodate stochastic volatility.}
\item{\code{\link{gFOBI}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{FOBI}} designed for time series with stochastic volatility.}
\item{\code{\link{gJADE}}}  {Function to solve the BSS problem. Algorithm is a generalization of \code{\link[JADE]{JADE}} designed for time series with stochastic volatility.}
\item{\code{\link{vSOBI}}}  {Function to solve the BSS problem. Algorithm is a variant of \code{\link[JADE]{SOBI}} algorithm and an alternative to \code{\link{FixNA}} to acommodate stochastic volatility.}}
Methods for ARMA models, like AMUSE and SOBI, and some non-stationary BSS methods for time series are implemented in the \code{\link[JADE]{JADE-package}}.
}

\author{
Markus Matilainen, Jari Miettinen, Klaus Nordhausen, Hannu Oja, Sara Taskinen

Maintainer: Markus Matilainen <markus.matilainen@utu.fi> 
}
\references{
Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New independent component analysis tools for time series}, Statistics & Probability Letters, 105, 80--87.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2016), \emph{ICA and Stochastic Volatility Models}, Proceedings of the XI International Conference on Computer Data Analysis and Modeling (CDAM'2016), 30--37.
}

\keyword{package}
\keyword{multivariate}
\keyword{ts}
