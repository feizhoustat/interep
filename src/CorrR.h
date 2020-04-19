#ifndef CorrR_h
#define CorrR_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::cube CorrR(arma::vec& y, arma::mat& x, arma::vec& beta, int n, arma::vec& k, char corre);
#endif
