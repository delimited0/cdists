// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
//#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;


// [[Rcpp::export]]
arma::cube cpropose_gwish(int n, arma::mat Sigma, double nu, arma::mat a, int npairs) {
//    std::default_random_engine generator;
    int d = Sigma.n_rows;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    double term1;
    double term2;
    double outterm1;
    double premult;
    double term4;
    double term3;
    double innerterm;
    double outerterm;
    double innerterm2;
    double outerterm2;
    
    
    int m;
    int p;
    
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
            
            if(a(m,p) != 0) {
                A(p,m) = as<double>(wrap(Rcpp::qnorm(Rcpp::runif(1), 0.0, 1.0, 1, 0))); 
            } else {
                term1 = 0.0;
                term2 = 0.0;
                outterm1 = 0.0;
                
                if(m == 0) {
                    premult = U(m,m) * A(m,m);
                    term4 = 0.0;
                    for(auto i = m; i < p; i++) {
                        term4 += U(i,p) * A(i,m);
                    }
                    A(p,m) = ( - term4 ) / U(p,p);
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

                    A(p,m) = (-term1 - term2 - term3) / (U(p,p) * U(m,m) * A(m,m));
                }
            }
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    return out;
}

