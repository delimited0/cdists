\name{propose_gwish}
\alias{propose_gwish}
\title{
Make proposals for the G-Wishart distribution
}
\usage{
propose_gwish(n, Psi, nu, adj)
}
\arguments{
  \item{n}{
Integer sample size.
}
  \item{Psi}{
A positive-definite matrix parameterizing the Wishart distribution.
}
  \item{nu}{
Degrees of freedom parameterizing the Wishart distribution.
}
  \item{adj}{
Symmetric matrix $X$ with 1 in $X_{ij}$ is there is an edge between nodes $i$ and $j$. The values along the diagonal do not affect the result.
}
}
\details{
Note that these draws do not correspond exactly to the G-Wishart distribution, in that a G-Wishart is like an inverse Wishart, while this function produced Wishart-like matrices.
}
\value{
An array of simulated matrices with zeros in the specified locations.
}
\author{\email{hbk5086@psu.edu}}

\examples{
    set.seed(2)
    d <- 5
    
    # Simulate adjacency matrix
    adj <- BDgraph::graph.sim(p = d, graph = "random", prob = .5)
    nu <- d + 10
    Psi <- LaplacesDemon::rinvwishart(nu = nu, S = diag(d) * 2)
    
    propose_gwish(10, Psi = Psi, nu = nu, adj = adj)    
}