
# compute R given the an initail estimate
Rf<-function(x){
  gap0<-sort(diff(c(0,sort(abs(x)))),decreasing = TRUE)
  return(gap0[1]/gap0[2])
}

# (1.3) AIC
LM_AIC<-function(y,x,beta,intercept=0){
  n<-length(y)
  SSE<-sum((y-cbind(1,x)%*%c(intercept,beta))^2)
  a<-n*log(SSE/n)+2*sum(beta!=0)
  return(a)
}

# (1.4) BIC
LM_BIC<-function(y,x,beta,intercept=0){
  n<-length(y)
  SSE<-sum((y-cbind(1,x)%*%c(intercept,beta))^2)
  b<-n*log(SSE/n)+log(n)*sum(beta!=0)
  return(b)
}

# (1.5) EBIC
LM_EBIC<-function(y,x,beta,gamma,intercept=0){
  n<-length(y)
  p<-dim(x)[2]
  SSE<-sum((y-cbind(1,x)%*%c(intercept,beta))^2)
  b<-n*log(SSE/n)+log(n)*sum(beta!=0)
  e<-b+2*gamma*log(choose(p,sum(beta!=0)))
  return(e)
}

LM_aEBIC<-function(y,x,beta,gamma,intercept=0){
  n<-length(y)
  p<-dim(x)[2]
  SSE<-sum((y-cbind(1,x)%*%c(intercept,beta))^2)
  b<-n*log(SSE/n)+log(n)*sum(beta!=0)
  e<-b+2*gamma*log(p)*sum(beta!=0)
  return(e)
}

# (1.6) Generalized cross vadiation
LM_GCV<-function(y,x,beta,intercept=0){
  n<-length(y)
  SSE<-sum((y-cbind(1,x)%*%c(intercept,beta))^2)
  g<-SSE/(n*(1-sum(beta!=0)/n)^2)
  return(g)
}
