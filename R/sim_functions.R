rcwish <- function(n, Psi, nu, a, b) {
    d <- nrow(Psi)
    npairs <- choose(d, 2)
    crcwish(n, Psi, nu, a, b, npairs)
}

propose_gwish <- function(n, Psi, nu, adj) {
    d <- nrow(Psi)
    npairs <- choose(d, 2)
    cpropose_gwish(n, Psi, nu, adj, npairs)
}

rtmvn <- function(n, Sigma, mu, a, b) {
    samp <- crtmvn(n = n, Sigma = Sigma, mu = mu, a = a, b = b)
    sweep(samp, 2, mu, "+")
}