#' @useDynLib interep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' fit generalized estimaitng equations with given tuning parameters
#'
#' This function makes predictions for generalized estimating equation with a given value of lambda.
#' Typical usage is to have the cv.interep function compute the optimal lambda, then provide it to
#' the interep function.
#' @importFrom stats gaussian
#' @importFrom MASS ginv
#' @param e matrix of environment factors.
#' @param g matrix of omics factors. In the case study, the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param beta0 the inital coefficient vector.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "a" as AR-1", "i" as "independence" and "e" as "exchangeable".
#' @param pmethod the penalization method. "mixed" refers to MCP penalty to individual main effects and group MCP penalty to interactions; "individual" means MCP penalty to all effects.
#' @param lam1 the tuning parameter lambda1 for individual predictors.
#' @param lam2 the tuning parameter lambda2 for interactions.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm.  The default value is 30
#' @return
#' \item{coef}{the coefficient vector.}
#' @references
#' Zhou, F., Ren, J., Li, G., Jiang, Y., Li, X., Wang, W.and Wu, C. (2019). Penalized variable selection for Lipid--environment interactions in the longitudinal lipidomics study.
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang,Y. and Wu, C. (2019). Semi-parametric Bayesian variable selection for Gene-Environment interactions.
#'
#' Wu, C., Zhou, F., Ren, J., Li, X., Jiang, Y., Ma, S. (2019). A Selective Review of Multi-Level Omics Data Integration Using Variable Selection.
#' \href{https://doi.org/10.3390/ht8010004}{\emph{High-Throughput}, 8(1)}
#'
#' Wu, C., Zhong, P.-S., and Cui, Y. (2018). Additive varying-coefficient model for nonlinear gene-environment interactions.
#' \href{https://doi.org/10.1515/sagmb-2017-0008}{\emph{ Statistical Applications in Genetics and Molecular Biology}, 17(2)}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y., Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' \href{https://doi.org/10.1002/sim.7518}{\emph{Statistics in Medicine}, 37:437–456}
#'
#' Jiang, Y., Huang, Y., Du, Y., Zhao, Y., Ren, J., Ma, S., & Wu, C. (2017). Identification of prognostic genes and pathways in lung adenocarcinoma using a Bayesian approach.
#' \href{https://www.researchgate.net/profile/Cen_Wu/publication/310828910_Identification_of_Prognostic_Genes_and_Pathways_in_Lung_Adenocarcinoma_Using_a_Bayesian_Approach/links/5cfbdac692851c874c5947f6/Identification-of-Prognostic-Genes-and-Pathways-in-Lung-Adenocarcinoma-Using-a-Bayesian-Approach.pdf}{\emph{Cancer Inform}, 1(7)}
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' \href{https://doi.org/10.1093/bib/bbu046}{\emph{Briefings in Bioinformatics}, 16(5), 873–883}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' \href{https://doi.org/10.1002/sim.6609}{\emph{Statistics in Medicine}, 34 (30): 4016–4030}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' \href{https://doi.org/10.1002/sim.6287}{\emph{Statistics in Medicine}, 33(28), 4988–4998}
#'
#' Wu, C. and Cui, Y. (2013). A novel method for identifying nonlinear gene–environment interactions in case–control association studies.
#' \href{https://doi.org/10.1007/s00439-013-1350-z}{\emph{Human Genetics}, 132(12):1413–1425}
#'
#' Wu, C. and Cui, Y. (2013). Boosting signals in gene–based association studies via efficient SNP selection.
#' \href{https://doi.org/10.1093/bib/bbs087}{\emph{Briefings in Bioinformatics}, 15(2):279–291}
#'
#' Wu, C., Li, S., and Cui, Y. (2012). Genetic Association Studies: An Information Content Perspective.
#' \href{https://doi.org/10.2174/138920212803251382}{\emph{Current Genomics}, 13(7),  566–573}
#'
#' @examples
#' data("dat")
#' e=dat$e
#' g=dat$z
#' y=dat$y
#' beta0=dat$coef
#' index=dat$index
#' b = interep(e, g, y,beta0,corre="e",pmethod="mixed",lam1=dat$lam1, lam2=dat$lam2,maxits=30)
#' b[abs(b)<0.05]=0
#' pos = which(b != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

interep <- function(e,g,y,beta0,corre,pmethod,lam1,lam2,maxits){
  q=dim(e)[2]
  n=dim(y)[1]
  k=dim(y)[2]
  p1=dim(g)[2]

  x=cbind(e,g)
  for (i in 1:p1) {
    for (j in 1:q) {
      x=cbind(x,e[,j]*g[,i])
    }
  }

  x=scale(x)

  #==========================================reformat the data===============================#
  data=reformat(k,y,x)
  y=data$y
  x=data$x

  p=dim(x)[2]
  k=rep(k,n)

  converge=F
  iter=0
  beta.new=beta0
  #=========================================PQIF======================================#
  while ((!converge) & (iter < maxits)) {
    beta = beta.new
    Score=ScoreU(n,k,y,x,p,beta,corre)
    U=Score$U
    dU=Score$dU

    E=penalty(x,n,p,q,beta,lam1,pmethod,p1,lam2)

    mat=dU + n*E
    mat[abs(mat)<0.000001]=0
    beta.new = beta + MASS::ginv(mat)%*%(U - n*E%*%beta)

    diff=mean(abs(beta.new-beta))
    converge = (diff < 1e-3)
    iter = iter+1
    cat("iter",iter,"diff",diff,"\n")
  }
  coef=beta.new
  return(coef)
}
