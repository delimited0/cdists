rcwish <- function(n, Psi, nu, a, b) {
    d <- nrow(Psi)
    npairs <- choose(d, 2)
    if(n == 1) {
        crcwish(n, Psi, nu, a, b, npairs)[,,1]
    } else {
        crcwish(n, Psi, nu, a, b, npairs)
    }
}

propose_gwish <- function(n, Psi, nu, adj) {
    d <- nrow(Psi)
    npairs <- choose(d, 2)
    if(n == 1) {
        cpropose_gwish(n, Psi, nu, adj, npairs)[,,1]
    } else {
        cpropose_gwish(n, Psi, nu, adj, npairs)   
    }
}

rtmvn <- function(n, Sigma, mu, a, b) {
    samp <- crtmvn(n = n, Sigma = Sigma, mu = mu, a = a, b = b)
    sweep(samp, 2, mu, "+")
}

rtmvt <- function(n, Sigma, mu, nu, a, b) {
    samp <- crtmvt(n = n, Sigma = Sigma, mu = mu, nu = nu, a = a, b = b)
    sweep(samp, 2, mu, "+")
}