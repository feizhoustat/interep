n=250
p=75
k=5
q=3
set.seed(123)
  y = matrix(rep(0,n*k),n,k)

  sig = matrix(0,p,p)

  for (i in 1: p)
  {
    for (j in 1: p)
    {
      sig[i,j] = 0.5^abs(i-j)
    }

  }
  x = MASS::mvrnorm(n,rep(0,p),sig)

  ###generate dummy variables
  dummy0 <- as.numeric(x[,1] <= -0.5)
  dummy0=dummy0+1
  dummy1 <- as.numeric(x[,1] > -0.5 & x[,1] <= 0)
  dummy1=dummy1*2+1
  dummy2 <- as.numeric(x[,1] > 0 & x[,1] <= 0.5)
  dummy2=dummy2*3+1

  x=cbind(dummy0,dummy1,dummy2,x)
  for (i in (q+1):(p+q)) {
    for (j in 1:q) {
      x=cbind(x,x[,j]*x[,i])
    }
  }

  x=scale(x)

  ll=0.4
  ul=0.8
  coef1=runif(q,ll,ul)
  #coef1=runif(q,-ul,-ll)
  coef2=runif(q,ll,ul)
  coef3=runif(q,ll,ul)
  coef4=runif(7,ll,ul)
  coef=c(coef4,coef1,coef2,coef3)
  mat=x[,c(1,2,3,5,7,10,15,(p+q+1):(p+q+3),(p+5*q+1):(p+5*q+3),(p+10*q+1):(p+10*q+3))]

  for(u in 1:k){
    y[,u] =  0.6+rowSums(coef*mat)
  }

  #AR(1) correlation matrix
  sig1 = matrix(0,k,k)
  diag(sig1)=1
  for (i in 1: k)
  {
    for (j in 1: k)
    {
      sig1[i,j] = 0.8^abs(i-j)
    }

  }
  error = MASS::mvrnorm(n,rep(0,k),sig1)

  y = y + error
  index=c(1,2,3,4,6,8,11,16,(p+q+2):(p+q+4),(p+5*q+2):(p+5*q+4),(p+10*q+2):(p+10*q+4))
  lam1=0.45
  lam2=1

  x1=cbind(data.frame(rep(1,n)),x)
  x1=data.matrix(x1)
  lasso.cv <- glmnet::cv.glmnet(x1,y[,3],alpha=1,nfolds=5)
  alpha <- lasso.cv$lambda.min/10  # lambda in the notes
  lasso.fit <- glmnet::glmnet(x1,y[,3],family="gaussian",alpha=1,nlambda=100)
  beta.new <- as.vector(stats::predict(lasso.fit, s=alpha, type="coefficients"))[-1]

  dat = list(y=y,e=x[,c(1:q)],z=x[,c((q+1):(p+q))],x=x,k=k,beta=beta.new,index=index,lam1=lam1,lam2=lam2)



