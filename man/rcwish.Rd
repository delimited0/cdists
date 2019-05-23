\name{rcwish}
\alias{rcwish}
\title{
Simulate from a Wishart distribution with constraints on the off-diagonals
}
\description{
The constraints on the off-diagonals $\rho$ take the form $a \leq \rho \leq b$ or $\rho = a$.
}
\usage{
rcwish(n, Psi, nu, a, b)
}
\arguments{
  \item{n}{integer sample size}
  \item{Psi}{A positive-definite matrix parameterizing the Wishart distribution.}
  \item{nu}{Degrees of freedom parameterizing the Wishart distribution.}
  \item{a}{Symmetrix matrix specifying the lower bounds on the off-diagonals. These can be set to -Inf. The values along the diagonal do not affect the result.}
  \item{b}{Symmetrix matrix specifying the upper bounds on the off-diagonals. These can be set to Inf. The values along the diagonal do not affect the result.}
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

a <- matrix(0, d, d)
b <- matrix(Inf, d, d)
rcwish(10, Psi = Psi, nu = nu, a = a, b = b)
}
