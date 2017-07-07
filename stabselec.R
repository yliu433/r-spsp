

# stability for lasso =========

stab_lasso <- function(x,y,nbootstrap=100,K=100,plotme=FALSE,thres=0.9)
{
  # Stability selection in the spirit of Meinshausen&Buhlman
  # JP Vert, 14/9/2010
  
  # x is the n*p design matrix, y the n*1 variable to predict
  require(glmnet)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  
  fl <- glmnet(x,y,family="gaussian",alpha=1,intercept=F,
         standardize = F,nlambda = K)
  
  lam <- fl$lambda
  freq <- matrix(0,length(lam),p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    # xs <- t(t(x)*runif(p,alpha,1))
    xs <- x
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    r <- glmnet(xs[i1,],y[i1],family="gaussian",alpha=1,intercept=F,
                standardize = F,lambda = lam)
    BETA <- t(as.matrix(r$beta))
    freq <- freq + abs(sign(BETA))
    
    r <- glmnet(xs[i2,],y[i2],family="gaussian",alpha=1,intercept=F,
                standardize = F,lambda = lam)
    BETA <- t(as.matrix(r$beta))
    freq <- freq + abs(sign(BETA))
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="lasso iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
  
  nz <- which(result > thres)
  z <- which(result <= thres)
  beta_stab<-rep(0,p)
  
  if(length(nz)>=1){
    Xc<-x[,nz]
    
    if(length(nz)<=n){
      betac<-lm(Y~Xc-1)$coefficients
    }else{
      betac<-solve(t(Xc)%*%Xc+0.001*diag(length(nz)))%*%t(Xc)%*%Y
    }
    
    beta_stab[nz]<-betac
  }  
  
  list(beta_stab=beta_stab,nonzero=nz,zero=z,freq_max=result)
  
}

# stability for adalasso =========

stab_adalasso <- function(x,y,nbootstrap=100,K=100,plotme=FALSE,thres=0.9)
{
  # Stability selection in the spirit of Meinshausen&Buhlman
  # JP Vert, 14/9/2010
  
  # x is the n*p design matrix, y the n*1 variable to predict
  require(glmnet)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  
  cv_ridge  <-  cv.glmnet(x, y, family='gaussian', alpha=0, standardize=F,intercept=F)
  w3  <-  1/abs(matrix(coef(cv_ridge, s=cv_ridge$lambda.min)[-1]))^1
  
  
  fl <- glmnet(x,y,family="gaussian",alpha=1,intercept=F,
                standardize = F,penalty.factor = w3,nlambda = K)
  
  lam <- fl$lambda
  freq <- matrix(0,length(lam),p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    # xs <- t(t(x)*runif(p,alpha,1))
    xs <- x
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    

    
    
    r <- glmnet(x[i1,],y[i1],family="gaussian",alpha=1,intercept=F,
                 standardize = F,penalty.factor = w3,lambda = lam)
    BETA <- t(as.matrix(r$beta))
    freq <- freq + abs(sign(BETA))
    

    
    
    r <- glmnet(x[i2,],y[i2],family="gaussian",alpha=1,intercept=F,
                standardize = F,penalty.factor = w3,lambda = lam)
    BETA <- t(as.matrix(r$beta))
    freq <- freq + abs(sign(BETA))
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="adalasso iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
  
  nz <- which(result > thres)
  z <- which(result <= thres)
  beta_stab<-rep(0,p)
  
  if(length(nz)>=1){
    Xc<-x[,nz]
    
    if(length(nz)<=n){
      betac<-lm(Y~Xc-1)$coefficients
    }else{
      betac<-solve(t(Xc)%*%Xc+0.001*diag(length(nz)))%*%t(Xc)%*%Y
    }
    
    beta_stab[nz]<-betac
  }  
  
  list(beta_stab=beta_stab,nonzero=nz,zero=z,freq_max=result)
  
}




# stability for scad =========

stab_scad <- function(x,y,nbootstrap=100,K=100,plotme=FALSE,thres=0.9)
{
  # Stability selection in the spirit of Meinshausen&Buhlman
  # JP Vert, 14/9/2010
  
  # x is the n*p design matrix, y the n*1 variable to predict
  require(ncvreg)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  
  fl <- ncvreg(x,y,family = "gaussian",penalty="SCAD",nlambda = K)
  
  lam <- fl$lambda
  freq <- matrix(0,length(lam),p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    # xs <- t(t(x)*runif(p,alpha,1))
    xs <- x
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    r <- ncvreg(x[i1,],y[i1],family = "gaussian",penalty="SCAD",lambda = lam)
    
    BETA <- t(as.matrix(r$beta[-1,]))
    freq <- freq + abs(sign(BETA))
    
    r <- ncvreg(x[i2,],y[i2],family = "gaussian",penalty="SCAD",lambda = lam)
    
    BETA <- t(as.matrix(r$beta[-1,]))
    freq <- freq + abs(sign(BETA))
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="scad iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
  
  nz <- which(result > thres)
  z <- which(result <= thres)
  beta_stab<-rep(0,p)
  
  if(length(nz)>=1){
    Xc<-x[,nz]
    
    if(length(nz)<=n){
      betac<-lm(Y~Xc-1)$coefficients
    }else{
      betac<-solve(t(Xc)%*%Xc+0.001*diag(length(nz)))%*%t(Xc)%*%Y
    }
    
    beta_stab[nz]<-betac
  }  
  
  list(beta_stab=beta_stab,nonzero=nz,zero=z,freq_max=result)
  
}

# stability for mcp =========

stab_mcp <- function(x,y,nbootstrap=100,K=100,plotme=FALSE,thres=0.9)
{
  # Stability selection in the spirit of Meinshausen&Buhlman
  # JP Vert, 14/9/2010
  
  # x is the n*p design matrix, y the n*1 variable to predict
  require(ncvreg)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  
  fl <- ncvreg(x,y,family = "gaussian",penalty="MCP",nlambda = K)
  
  lam <- fl$lambda
  freq <- matrix(0,length(lam),p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    # xs <- t(t(x)*runif(p,alpha,1))
    xs <- x
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    r <- ncvreg(x[i1,],y[i1],family = "gaussian",penalty="MCP",lambda = lam)
    
    BETA <- t(as.matrix(r$beta[-1,]))
    freq <- freq + abs(sign(BETA))
    
    r <- ncvreg(x[i2,],y[i2],family = "gaussian",penalty="MCP",lambda = lam)
    
    BETA <- t(as.matrix(r$beta[-1,]))
    freq <- freq + abs(sign(BETA))
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="mcp iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
  
  nz <- which(result > thres)
  z <- which(result <= thres)
  beta_stab<-rep(0,p)
  
  if(length(nz)>=1){
    Xc<-x[,nz]
    
    if(length(nz)<=n){
      betac<-lm(Y~Xc-1)$coefficients
    }else{
      betac<-solve(t(Xc)%*%Xc+0.001*diag(length(nz)))%*%t(Xc)%*%Y
    }
    
    beta_stab[nz]<-betac
  }  
  
  list(beta_stab=beta_stab,nonzero=nz,zero=z,freq_max=result)
  
}

