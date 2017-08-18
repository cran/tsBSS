\name{plot.tssdr}
\alias{plot.tssdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting an object of class tssdr}
\description{
Plots the response and the estimated directions (sources) resulting from a tssdr method. The multivariate time series is a time series object, which is passed on to \code{plot.ts}.
}
\usage{
\method{plot}{tssdr}(x, main = "The response and the directions", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class tssdr.}
  \item{main}{A title for the time series plot.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\author{Markus Matilainen}
\seealso{
\code{\link{plot.ts}}, \code{\link[tsBSS]{tssdr}}
}
\keyword{methods}
