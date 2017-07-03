\name{print.tssdr}
\alias{print.tssdr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing an object of class tssdr}
\description{
Prints an object of class tssdr. It prints the signal separation matrix W, matrix L to decide which lags and directions are important, and vector or scalar k, i.e. the lag(s) used.
}
\usage{
\method{print}{tssdr}(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class tssdr.}
  \item{digits}{Chosen number of digits.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}
\author{Markus Matilainen}
\seealso{
\code{\link[tsBSS]{tssdr}}
}
\keyword{methods}
