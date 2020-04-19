#ifndef reformat_h
#define reformat_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List reformat(arma::mat& y, arma::mat& x);

#endif
