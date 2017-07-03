\name{plot.summary.tssdr}
\alias{plot.summary.tssdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting an Object of class summary.tssdr}
\description{
Plots the response and the estimated directions (sources) resulting from a summary of a tssdr method. The directions are the chosen directions from tssdr method. The multivariate time series is a time series object, which is passed on to \code{plot.ts}. 
}
\usage{
\method{plot}{summary.tssdr}(x, main = "The response and the chosen directions", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class summary.tssdr.}
  \item{main}{Title of the time series plot.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\author{Markus Matilainen}
\seealso{
\code{\link{plot.ts}}, \code{\link[tsBSS]{summary.tssdr}}
}
\keyword{methods}
