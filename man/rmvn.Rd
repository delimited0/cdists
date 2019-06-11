\name{rmvn}
\alias{rmvn}
\title{
Simulate multivariate normal data
}
\description{
This function makes ordinary multivariate normal draws, for comparison to the truncated counterpart.
}
\usage{
rmvn(n, Sigma, mu)
}
\arguments{
  \item{n}{
Integer sample size.
}
  \item{Sigma}{
Positive-definite covariance matrix parameterizing the multivariate normal distribution.
}
  \item{mu}{
Mean vector parameterizing the multivariate normal distribution.
}
}
\value{
An n by d matrix, where d is the dimension.
}
\author{\email{hbk5086@psu.edu}}
\examples{
set.seed(3)
d <- 5
Sigma <- diag(d) * 2
offdiags <- matrix(sample(c(.8, -.8), size = d * d, replace = TRUE), nrow = d, ncol = d)
Sigma[lower.tri(Sigma)] <- offdiags[lower.tri(offdiags)]
Sigma <- Sigma + t(Sigma) + diag(d) * 8

mu <- rep(0.5, d)

rmvn(100, Sigma, mu)

}