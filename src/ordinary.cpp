// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
//#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::export]]
arma::cube cr_ordinary_wish(int n, arma::mat Sigma, double nu, int npairs) {
    int d = Sigma.n_rows;
    int m;
    int p;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    
    arma::rowvec dp;
    arma::mat dimpairs(npairs, 2, arma::fill::none);
    int counter = 0;
    for(double k = 0; k < (d-1); k++) {
        for(double i = 0; i <= k; i++) {
            dp = { i, k + 1 };
            dimpairs.row(counter) = dp;
            counter += 1;
        }
    }
    
    for(auto iteration = 0; iteration < n; iteration++) {
        for(auto i = 0; i < d; i++) {
            A(i,i) = sqrt(as<double>(wrap(Rcpp::rchisq(1, nu - i))));
        }
            
        for(auto idx = 0; idx < npairs; idx++) {
            m = dimpairs(idx,0);
            p = dimpairs(idx,1);
            A(p,m) = as<double>(wrap(Rcpp::qnorm(Rcpp::runif(1), 0.0, 1.0, 1, 0)));
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    
    return out;
}

// [[Rcpp::export]]
arma::mat cr_ordinary_mvn(int n, arma::mat Sigma, arma::colvec mu) {
    int d = mu.size();
    arma::mat L = chol(Sigma, "lower");
    arma::mat z(d, n, arma::fill::zeros);
    
    for(auto iteration = 0; iteration < n; iteration++){
        for(auto i = 0; i < d; i++) {
            z(i, iteration) = as<double>(wrap(Rcpp::qnorm(Rcpp::runif(1), 0.0, 1.0, 1, 0)));
        }
    }
    return((L * z).t());
}
