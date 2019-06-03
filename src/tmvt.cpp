// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

using namespace boost::math;
using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
arma::mat crtmvt(int n, arma::mat Sigma, arma::colvec mu, double nu, arma::colvec a, arma::colvec b) {
    std::default_random_engine generator;
    int d = mu.size();
    arma::mat L = chol(Sigma, "lower");
    double A;
    arma::mat z(d, n, arma::fill::zeros);
    double lb;
    double ub;
    normal standard_normal;
    std::uniform_real_distribution<double> standard_uniform(0.0, 1.0);
    
    chi_squared chidist(nu);
    A = sqrt(nu / quantile(chidist, standard_uniform(generator)));
    arma::mat AL = A * L;

    for(auto iteration = 0; iteration < n; iteration++){
        for(auto i = 0; i < d; i++) {
            lb = a(i) - mu(i);
            ub = b(i) - mu(i);

            if(i == 0) {
                lb /= AL(0,0);
                ub /= AL(0,0);
            } else {
                for(auto j = 0; j < i; j++) {
                    lb -= (AL(i,j) * z(j, iteration));
                    ub -= (AL(i,j) * z(j, iteration));
                }
                lb /= AL(i,i);
                ub /= AL(i,i);
            }
            std::uniform_real_distribution<double> runif(cdf(standard_normal, lb), cdf(standard_normal, ub));
            z(i, iteration) = quantile(standard_normal, runif(generator));
        }
    }
    return((AL * z).t());
}

