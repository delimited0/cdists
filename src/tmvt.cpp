// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
arma::mat crtmvt(int n, arma::mat Sigma, arma::colvec mu, double nu, arma::colvec a, arma::colvec b, arma::colvec scale_term) {
    int d = mu.size();
    int use_lower;
    arma::mat L = chol(Sigma, "lower");
    arma::mat z(d, n, arma::fill::zeros);
    double lb, log_lbnew;
    double ub, log_ubnew;
    
    double max_const;
    NumericVector u(1);
    NumericVector ru(1);
    
    arma::mat AL;
    arma::mat out(d, n, arma::fill::none);
    
    for(auto iteration = 0; iteration < n; iteration++){
        AL = scale_term(iteration) * L;
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
        out.col(iteration) = AL * (z.col(iteration));
    }
    return(out.t());
}

