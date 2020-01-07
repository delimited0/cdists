// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
//#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::export]]
arma::mat crtmvn(int n, arma::mat Sigma, arma::colvec mu, arma::colvec a, arma::colvec b) {
    int d = mu.size();
    int use_lower;
    arma::mat L = arma::chol(Sigma, "lower");
    arma::mat z(d, n, arma::fill::zeros);
    double lb, log_lbnew;
    double ub, log_ubnew;
    double max_const;
    NumericVector u(1);
    NumericVector ru(1);
    
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
            if(ub > 0){
                use_lower = 0;
            }
            else{
                use_lower = 1;
            }
      		log_lbnew = R::pnorm(lb, 0.0, 1.0, use_lower, 1);
            log_ubnew = R::pnorm(ub, 0.0, 1.0, use_lower, 1);
            max_const = max(NumericVector::create(log_lbnew, log_ubnew));
            ru = Rcpp::runif(1);
            u = log(exp((log(1 - ru) + log_lbnew  - max_const)) + 
                        exp((log(ru) + log_ubnew - max_const))) + max_const;
            z(i, iteration) = as<double>(wrap(Rcpp::qnorm(u, 0.0, 1.0, use_lower, 1))); 
        }
    }
    return((L * z).t());
}
