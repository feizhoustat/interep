
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include"dmcp.h"
#include<Rmath.h>
#include<iostream>
#include<stdio.h>

using namespace Rcpp;
using namespace arma;
using namespace R;

//' This function obtains the first derivative function of MCP (Minimax Concave Penalty)
//' @param theta a coefficient vector.
//' @param lambda the tuning parameter.
//' @param gamma the regularization parameter in MCP (Minimax Concave Penalty).
//' It balances between the unbiasedness and concavity of MCP.
//' @return the first derivative of MCP function.
//' @details
//' Rigorously speaking, the regularization parametre \eqn{\gamma} needs to be obtained via a data-driven approach.
//' Published studies suggest experimenting with a few values, such as 1.8, 3, 4.5, 6, and 10, then fixing its value. In our numerical
//' study, we have examined this sequence and found that the results are not sensitive to the choice of value of \eqn{\gamma},
//' and set the value at 3. In practice, to be prudent, values other than 3 should also be investigated.
//'
//' @examples
//' theta=runif(20,-5,5)
//' lambda=1
//' gamma=3
//' dmcp(theta,lambda,gamma)
//' @export
// [[Rcpp::export()]]
arma::vec dmcp(arma::vec& theta, double lambda, double gamma){
  int p=theta.size();
  arma::vec b(p,arma::fill::zeros),
      b1(p,arma::fill::zeros);
  for(int i=0; i<p; i++){
    if(std::abs(theta(i))<=lambda*gamma){
      b1(i)=1;
    }
    b(i)=(lambda-std::abs(theta(i))/gamma)*b1(i);
  }

  return b;
}


