% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{summary.fut}
\alias{summary.fut}
\title{Detailed summary of the futility design}
\usage{
\method{summary}{fut}(object, ...)
}
\arguments{
\item{object}{An object returned by \code{\link{fut}}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class \code{summary.fut} with components:
\item{t}{A \eqn{K}-dimensional vector of information times.}
\item{b}{A \eqn{K}-dimensional vector of B-values at \code{t}.}
\item{z}{A \eqn{K}-dimensional vector of z-values at \code{t}.}
\item{type2}{A \eqn{K}-dimensional vector of beta spent at \code{t}.}
\item{loss}{A \eqn{K}-dimensional vector of power loss at \code{t}.}
\item{ess}{Expected sample size at \eqn{H_0}.}
\item{...}{}
}
\description{
Provide key information about the futility design, including B-/z-values, beta (type II error)
spent, and power loss at each futility look as well the the sample size distribution under
the null hypothesis.
}
\examples{
# see example for fut
}
\seealso{
\code{\link{fut}}, \code{\link{print.fut}}, \code{\link{print.summary.fut}}.
}
\keyword{fut}
