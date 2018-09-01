#' This function obtains the first derivative function of MCP (Minimax Concave Penalty)
#' @param theta a coefficient vector.
#' @param lambda the tuning parameter.
#' @param gamma the regularization parameter in MCP (Minimax Concave Penalty).
#' It balances between the unbiasedness and concavity of MCP.
#' @return the first derivative of MCP function.
#' @details
#' Rigorously speaking, the regularization parametre \eqn{\gamma} needs to be obtained via a data-driven approach.
#' Published studies suggest experimenting with a few values, such as 1.8, 3, 4.5, 6, and 10, then fixing its value. In our numerical
#' study, we have examined this sequence and found that the results are not sensitive to the choice of value of \eqn{\gamma},
#' and set the value at 3. In practice, to be prudent, values other than 3 should also be investigated.
#'
#' @examples
#' theta=runif(20,-5,5)
#' lambda=1
#' dmcp(theta,lambda,gamma=3)
#' @export

dmcp <- function(theta,lambda,gamma=3){
  #length of parameter
  p<-length(theta)

  #create vector of zeros
  b1<-rep(0,p)
  #if theta is less than gamma*lambda set it to 1
  b1[abs(theta)<(gamma*lambda)]<-1
  (lambda-theta/gamma)*b1
}
