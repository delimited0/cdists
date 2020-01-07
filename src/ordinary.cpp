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

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
arma::cube cr_ordinary_wish(int n, arma::mat Sigma, double nu, int npairs) {
    std::default_random_engine generator{(unsigned) randWrapper(INT_MAX)};
    int d = Sigma.n_rows;
    int m;
    int p;
    arma::mat U = chol(Sigma, "upper");
    arma::mat A(d, d, arma::fill::zeros);
    arma::cube out(d, d, n, arma::fill::zeros);
    
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
            A(p,m) = quantile(standard_normal, standard_uniform(generator));
        }
        out.slice(iteration) = (U.t() * A * A.t() * U);
    }
    
    return out;
}

// [[Rcpp::export]]
arma::mat cr_ordinary_mvn(int n, arma::mat Sigma, arma::colvec mu) {
    std::default_random_engine generator;
    int d = mu.size();
    arma::mat L = chol(Sigma, "lower");
    arma::mat z(d, n, arma::fill::zeros);
    
    normal standard_normal;
    std::uniform_real_distribution<double> standard_uniform(0.0, 1.0);

    for(auto iteration = 0; iteration < n; iteration++){
        for(auto i = 0; i < d; i++) {
            z(i, iteration) = quantile(standard_normal, standard_uniform(generator));
        }
    }
    return((L * z).t());
}

