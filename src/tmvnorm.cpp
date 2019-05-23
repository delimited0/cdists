// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <random>
#include <boost/math/distributions/normal.hpp>

using namespace boost::math;
using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
arma::mat crtmvn(int n, arma::mat Sigma, arma::colvec mu, arma::colvec a, arma::colvec b) {
    std::default_random_engine generator;
    int d = mu.size();
    arma::mat L = chol(Sigma, "lower");
    arma::mat z(d, n, arma::fill::zeros);
    double lb;
    double ub;
    normal standard_normal;

    for(auto iteration = 0; iteration < n; iteration++){
        for(auto i = 0; i < d; i++) {
            lb = a(i) - mu(i);
            ub = b(i) - mu(i);

            if(i == 0) {
                lb /= L(0,0);
                ub /= L(0,0);
            } else {
                for(auto j = 0; j < i; j++) {
                    lb -= (L(i,j) * z(j, iteration));
                    ub -= (L(i,j) * z(j, iteration));
                }
                lb /= L(i,i);
                ub /= L(i,i);
            }
            std::uniform_real_distribution<double> runif(cdf(standard_normal, lb), cdf(standard_normal, ub));
            z(i, iteration) = quantile(standard_normal, runif(generator));
        }
    }
    return((L * z).t());
}

