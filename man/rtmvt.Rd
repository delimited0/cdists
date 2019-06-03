\name{rtmvt}
\alias{rtmvt}
\title{
Simulate from the truncated multivariate Student's t distribution
}
\description{
This algorithm makes exact draws from the truncated multivariate Student's t distribution.
}
\usage{
rtmvt(n, Sigma, mu, nu, a, b)
}
\arguments{
  \item{n}{
Integer sample size.
}
  \item{Sigma}{
Positive-definite covariance matrix parameterizing the multivariate Student's t distribution.
}
  \item{mu}{
Non-centrality vector parameterizing the multivariate Student's t distribution.
}
\item{nu}{
    Degrees of freedom.
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
Sigma[lower.tri(Psi)] <- offdiags[lower.tri(offdiags)]
Sigma <- Sigma + t(Sigma) + diag(d) * 8

a <- rep(0, d)
b <- rep(Inf, d)
mu <- rep(0.5, d)
nu <- 50

rtmvt(100, Sigma, mu, nu, a, b)

}