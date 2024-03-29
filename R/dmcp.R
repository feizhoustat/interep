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
#' and set the value at 3. In practice, to be prudent, values other than 3 should also be investigated. Similar discussions can be found
#' in the references below.
#' @references
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang, Y. and Wu, C. (2019). Robust network-based regularization and variable selection for high-dimensional genomic data in cancer prognosis.
#' {\emph{Genetic epidemiology}, 43(3), 276-291} \doi{10.1002/gepi.22194}
#'
#' Ren, J., Jung, L., Du, Y., Wu, C., Jiang, Y. and Liu, J. (2019). regnet: Network-Based Regularization for Generalized Linear Models.
#' \href{https://cran.r-project.org/package=regnet}{\emph{R package}, version 0.4.0}
#'
#' Wu, C., Zhang, Q., Jiang, Y. and Ma, S. (2018). Robust network-based analysis of the associations between (epi) genetic measurements.
#' {\emph{Journal of multivariate analysis}, 168, 119-130} \doi{10.1016/j.jmva.2018.06.009}
#'
#' Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y. and Wu, C. (2017). Network-based regularization for high dimensional SNP data in the case–control study of Type 2 diabetes.
#' {\emph{BMC genetics}, 18(1), 44} \doi{10.1186/s12863-017-0495-5}
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
  b1[abs(theta)<=(gamma*lambda)]<-1
  (lambda-theta/gamma)*b1
}
