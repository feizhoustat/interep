#ifndef dmcp_h
#define dmcp_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec dmcp(arma::vec& theta, double lambda, double gamma);

#endif
