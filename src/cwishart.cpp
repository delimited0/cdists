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
arma::cube crcwish(int n, arma::mat Sigma, double nu, arma::mat a, arma::mat b, int npairs) {
    std::default_random_engine generator;
    int d = Sigma.n_rows;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    double lb;
    double ub;
    double term1;
    double term2;
    double outterm1;
    double premult;
    double term4;
    double outterm4;
    double denom;
    //double val;
    int m;
    int p;
    normal standard_normal;
    std::uniform_real_distribution<double> standard_uniform(0.0, 1.0);
    
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
            chi_squared chidist(nu-i);
            A(i,i) = sqrt(quantile(chidist, standard_uniform(generator)));
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
                outterm4 = premult * term4;
                denom = premult * U(p,p);
                
                if(a(m,p) == b(m,p)) {
                    A(p,m) = (a(m,p) - outterm4) / denom;
                } else {
                    lb = (a(m,p) - outterm4) / denom;
                    ub = (b(m,p) - outterm4) / denom;
                    std::uniform_real_distribution<double> runif(cdf(standard_normal, lb), cdf(standard_normal, ub));
                    A(p,m) = quantile(standard_normal, runif(generator));
                }
            } else {
                for(auto k = 0; k < m; k++) {
                    for(auto j = k; j <= m; j++) {
                        term1 += U(j,m) * A(j,k);
                    }
                    for(auto i = k; i <= p; i++) {
                        term2 += U(i,p) * A(i,k);
                    }
                    outterm1 += (term1 * term2);
                }
                
                premult = U(m,m) * A(m,m);
                
                term4 = 0.0;
                for(auto i = m; i < p; i++) {
                    term4 += U(i,p) * A(i,m);
                }
                outterm4 = premult * term4;
                denom = premult * U(p,p);
                
                if(a(m,p) == b(m,p)) {
                    A(p,m) = (a(m,p) - outterm1 - outterm4) / denom;
                } else {
                    lb = (a(m,p) - outterm1 - outterm4) / denom;
                    ub = (b(m,p) - outterm1 - outterm4) / denom;
                    std::uniform_real_distribution<double> runif(cdf(standard_normal, lb), cdf(standard_normal, ub));
                    A(p,m) = quantile(standard_normal, runif(generator));
                }
            }
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    
    return out;
}

// [[Rcpp::export]]
arma::cube cpropose_gwish(int n, arma::mat Sigma, double nu, arma::mat a, int npairs) {
    std::default_random_engine generator;
    int d = Sigma.n_rows;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    double term1;
    double term2;
    double outterm1;
    double premult;
    double term4;
    double outterm4;
    double denom;
    int m;
    int p;
    normal standard_normal;
    std::uniform_real_distribution<double> standard_uniform(0.0, 1.0);
    
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
            chi_squared chidist(nu-i);
            A(i,i) = sqrt(quantile(chidist, standard_uniform(generator)));
        }
        
        for(auto idx = 0; idx < npairs; idx++) {
            m = dimpairs(idx,0);
            p = dimpairs(idx,1);
            
            if(a(m,p) != 0) {
                std::uniform_real_distribution<double> runif(0.0, 1.0);
                A(p,m) = quantile(standard_normal, runif(generator));
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
                    outterm4 = premult * term4;
                    denom = premult * U(p,p);
                    A(p,m) = (-outterm4) / denom;
                } else {
                    for(auto k = 0; k < m; k++) {
                        for(auto j = k; j <= m; j++) {
                            term1 += U(j,m) * A(j,k);
                        }
                        for(auto i = k; i <= p; i++) {
                            term2 += U(i,p) * A(i,k);
                        }
                        outterm1 += (term1 * term2);
                    }
                    
                    premult = U(m,m) * A(m,m);
                    
                    term4 = 0.0;
                    for(auto i = m; i < p; i++) {
                        term4 += U(i,p) * A(i,m);
                    }
                    outterm4 = premult * term4;
                    denom = premult * U(p,p);
                    
                    A(p,m) = (-outterm1 - outterm4) / denom;
                }
            }
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    
    return out;
}

