\name{rwish}
\alias{rwish}
\title{
Simulate from a Wishart distribution
}
\description{
This function makes ordinary Wishart draws, for comparison to the constrained counterpart.
}
\usage{
rwish(n, Psi, nu)
}
\arguments{
  \item{n}{integer sample size}
  \item{Psi}{A positive-definite matrix parameterizing the Wishart distribution.}
  \item{nu}{Degrees of freedom parameterizing the Wishart distribution.}
}
\value{
An array of simulated matrices with the specified constraints.
}
\author{\email{hbk5086@psu.edu}}
\examples{
set.seed(1)
d <- 10
nu <- d + 5
Psi <- LaplacesDemon::rinvwishart(nu = nu, S = diag(d) * 2)

rwish(10, Psi = Psi, nu = nu)
}
