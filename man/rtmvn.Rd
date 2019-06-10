\name{rtmvn}
\alias{rtmvn}
\title{
Simulate from the truncated multivariate normal distribution
}
\description{
This algorithm makes exact draws from the truncated multivariate normal distribution.
}
\usage{
rtmvn(n, Sigma, mu, a, b)
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
  \item{a}{
Vector of lower bounds on the random draws in each dimension. These can be -Inf if desired.
}
  \item{b}{
Vector of upper bounds on the random draws in each dimension. These can be Inf if desired.
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

a <- rep(0, d)
b <- rep(Inf, d)
mu <- rep(0.5, d)

rtmvn(100, Sigma, mu, a, b)

}