// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
//#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;

// [[Rcpp::export]]
arma::cube crcwish(int n, arma::mat Sigma, double nu, arma::mat a, arma::mat b, int npairs) {
    int d = Sigma.n_rows;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    NumericVector u(1);
    NumericVector ru(1);
    double lb, log_lbnew;
    double ub, log_ubnew;
    double max_const;
    double outterm1;
    double premult;
    double term4;
    double term1;
    double term2;
    double term3;
    double innerterm;
    double outerterm;
    double innerterm2;
    double outerterm2;


    int m;
    int p;
    int use_lower;
    
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
            

            term1 = 0.0;
            term2 = 0.0;
            outterm1 = 0.0;
            
            if(m == 0) {
                premult = U(m,m) * A(m,m);
                term4 = 0.0;
                for(auto i = m; i < p; i++) {
                    term4 += U(i,p) * A(i,m);
                }
                
                if(a(m,p) == b(m,p)) {
                    A(p,m) = ( ( a(m,p) / premult ) - term4 ) / U(p,p);
                } else {
                    lb = ( ( a(m,p) / premult ) - term4 ) / U(p,p);
                    ub = ( ( b(m,p) / premult ) - term4 ) / U(p,p);
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
                    A(p,m) = as<double>(wrap(Rcpp::qnorm(u, 0.0, 1.0, use_lower, 1))); 
                }
            } else {
                term1 = 0.0;
                for(auto k = 0; k < m; k++) {
                    outerterm = 0.0;
                    for(auto j=0; j <=k; j++) {
                        innerterm = 0.0;
                        for(auto i = j; i <= p; i++) {
                            innerterm += (A(i,j) * U(i,p));
                        }
                        outerterm += (innerterm * A(k,j));
                    }
                    term1 += U(k,m) * outerterm;
                }

                outerterm2 = 0.0;
                for(auto j=0; j < m; j++) {
                    innerterm2 = 0.0;
                    for(auto i=j; i <= p; i++) {
                        innerterm2 += (A(i,j) * U(i,p));
                    }
                    outerterm2 += (innerterm2 * A(m,j));
                }
                term2 = U(m,m) * outerterm2;

                term3 = 0.0;
                for(auto i = m; i < p; i++) {
                    term3 += (A(i,m) * U(i,p));
                }
                term3 *= (A(m,m) * U(m,m));
                
                if(a(m,p) == b(m,p)) {
                    A(p,m) = (a(m,p) - term1 - term2 - term3) / (U(p,p) * U(m,m) * A(m,m));
                } else {
                    lb = (a(m,p) - term1 - term2 - term3) / (U(p,p) * U(m,m) * A(m,m));
                    ub = (b(m,p) - term1 - term2 - term3) / (U(p,p) * U(m,m) * A(m,m));
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
                    A(p,m) = as<double>(wrap(Rcpp::qnorm(u, 0.0, 1.0, use_lower, 1))); 
                }
            }
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    
    return out;
}
