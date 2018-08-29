#' This function obtains the correlation matrix of longitudinal data
#' @importFrom stats family gaussian sd uniroot
#' @param x a matrix of predictors, consisting of omics and environment factors, as well as their interactions. In the case study,
#'  the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param beta the coefficient vector.
#' @param n the sample size.
#' @param k the number of repeated measurement.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "AR-1", "independece" and "exchangeable".
#' @return the correlation matrix.
#' @export

CorrR <- function(y,x,beta,n,k,corre){
  aindex=cumsum(k)
  index=c(0,aindex[-length(aindex)])
  mu=x%*%beta
  res=(as.vector(y)-mu)/sd(y)

  if(corre == "independence"){
    alfa_hat=0
  }
  else if(corre == "exchangeable"){

    sum5<-0
    sum6<-0
    for ( i in  1:n)           {
      for ( j in  1:k[i])       {
        for ( jj in 1:k[i])       {
          if( j>jj)  {
            if(abs(j-jj)==1){
              #cat("i",i,"j",j,"jj",jj,"\n")
              sum7<-res[j+index[i]]*res[jj+index[i]]
              sum5<-sum5+sum7
              #cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")
            }

          }
        }
      }
      sum8<-(k[i]-1)
      sum6<-sum6+sum8
    } #i
    alfa_hat<-sum5/sum6

  }

  else if (corre == "AR-1") {
    sum5<-0
    sum6<-0
    for ( i in  1:n)           {
      for ( j in  1:k[i])       {
        for ( jj in 1:k[i])       {
          if( j>jj)  {

              #cat("i",i,"j",j,"jj",jj,"\n")
              sum7<-res[j+index[i]]*res[jj+index[i]]
              sum5<-sum5+sum7
              #cat("i",i,"j",j,"jj",jj,"sum7",sum7,"sum5", sum5, "\n")


          }
        }
      }
      sum8<-(k[i]-1)
      sum6<-sum6+sum8
    } #i
    alfa_hat<-uniroot(function(x) x^4+2*x^3+3*x^2+4*x-sum5/n, c(0,1), tol = 0.001)$root
  }


  maxclsz=max(k)

  Rhat<-array(0,c(maxclsz,maxclsz,n))

  for (i in 1:n) {
    cor1=matrix(0,k[i],k[i])
    if (corre == "independence") {
      cor1 = diag(k[i])
    }
    else if (corre == "exchangeable") {
      for (t1 in 1:k[i]) {
        for (t2 in 1:k[i]) {
          if (t1 != t2) {
            cor1[t1,t2] = alfa_hat
          }
          else {
            cor1[t1,t2]=1
          }
        }
      }
    }
    else if (corre == "AR-1"){
      for (t1 in 1:k[i]) {
        for (t2 in 1:k[i]) {
          cor1[t1,t2]<-alfa_hat^abs(t1-t2)
        }
      }
    }


    Rhat[1:k[i],1:k[i],i]<-cor1
  }

  return(Rhat)
}
