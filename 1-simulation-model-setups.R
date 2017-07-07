
require(MASS)
#### different Model settings

if(M==1){
  n <- 50;p <- 100;sigma <- sigma0
  beta <- rep(0,p);nonzero <- c(1,2,5);zero <- setdiff(1:p,nonzero)
  beta[nonzero] <- c(3,1.5,2)
  Sigma <- 0.5^(abs(outer(1:p,1:p,"-")))
  
  # restore the data
  Xdata <- matrix(0,n*N,p)  ### Coviates
  Ydata <- numeric(n*N)     ### Response
  for (j in 1:N){
    X <- mvrnorm(n,rep(0,p),Sigma)
    error <- rnorm(n,0,sigma)
    
    X <- apply(X,2,scale)*sqrt(n)/sqrt(n-1)
    error <- error-mean(error)
    
    Y <- X%*%beta+error
    Xdata[(1+(j-1)*n):(j*n),] <- X
    Ydata[(1+(j-1)*n):(j*n)] <- Y
  }
  
}else if(M==2){
  n <- 50;p <- 1000;sigma <- sigma0
  beta <- rep(0,p);nonzero <- c(1,2,5);zero <- setdiff(1:p,nonzero)
  beta[nonzero] <- c(3,1.5,2)
  Sigma <- 0.5^(abs(outer(1:p,1:p,"-")))
  
  # restore the data
  Xdata <- matrix(0,n*N,p)  ### Coviates
  Ydata <- numeric(n*N)     ### Response
  for (j in 1:N){
    X <- mvrnorm(n,rep(0,p),Sigma)
    error <- rnorm(n,0,sigma)
    
    X <- apply(X,2,scale)*sqrt(n)/sqrt(n-1)
    error <- error-mean(error)
    
    Y <- X%*%beta+error
    Xdata[(1+(j-1)*n):(j*n),] <- X
    Ydata[(1+(j-1)*n):(j*n)] <- Y
  }
  
}else if(M==3){
  n <- 50;p <- 100;sigma <- sigma0
  beta <- rep(0,p);nonzero <- 1:6;zero <- setdiff(1:p,nonzero)
  beta[nonzero] <- c(3,3,-2,3,3,-2)
  Sigma1 <- matrix(0.9,3,3);diag(Sigma1) <- 1
  #Sigma2 <- diag(p-10)
  
  # restore the data
  Xdata <- matrix(0,n*N,p)  ### Coviates
  Ydata <- numeric(n*N)     ### Response
  for (j in 1:N){
    X1 <- mvrnorm(n,rep(0,3),Sigma1)
    X2 <- mvrnorm(n,rep(0,3),Sigma1)
    X3 <- matrix(rnorm((p-6)*n),n,p-6)
    
    X <- cbind(X1,X2,X3)
    error <- rnorm(n,0,sigma)
    
    X <- apply(X,2,scale)*sqrt(n)/sqrt(n-1)
    error <- error-mean(error)
    
    Y <- X%*%beta+error
    Xdata[(1+(j-1)*n):(j*n),] <- X
    Ydata[(1+(j-1)*n):(j*n)] <- Y
  }
  
}else if(M==4){
  n <- 50;p <- 100;sigma <- min(sigma0/4,1)
  beta <- rep(0,p);nonzero <- 1:5;zero <- setdiff(1:p,nonzero)
  beta[nonzero] <- c(1,-1.25,0.75,-0.95,1.5)
  # restore the data
  Xdata <- matrix(0,n*N,p)  ### Coviates
  Ydata <- numeric(n*N)     ### Response
  for (j in 1:N){
    
    X <- matrix(rnorm(p*n),n,p)
    
    error <- rnorm(n,0,sigma)
    
    X <- apply(X,2,scale)*sqrt(n)/sqrt(n-1)
    error <- error-mean(error)
    
    Y <- X%*%beta+X[,1]*X[,2]+error
    Xdata[(1+(j-1)*n):(j*n),] <- X
    Ydata[(1+(j-1)*n):(j*n)] <- Y
  }
  
}