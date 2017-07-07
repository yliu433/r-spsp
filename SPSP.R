
### Select the informative variables by the SPSP algorithm

SPSP<-function(Y,X,BETA,init=1,R=NULL,standardize=FALSE){
  
  n<-dim(X)[1] # size
  p<-dim(X)[2] # dimension
  BETA<-as.matrix(BETA) 
  # BETA should be a p by k matrix;should be sparser and sparser
  # each column corrspends to a lambda, and lambda gets larger and larger
  K<-dim(BETA)[2] # the number of the tuning parameter lambda
  
  if(is.null(R)){
    #### select the default the tuning parameter if it is not given
    gap0<-sort(diff(c(0,sort(abs(BETA[,init])))),decreasing = TRUE)
    # sorted adjacent distances 
    R<-ifelse(gap0[2]!=0,as.numeric(gap0[1]/gap0[2]),as.numeric(gap0[1]))   
  }

  ### only consider the rows with at least one non-zero beta
  colsbeta<-colSums(BETA)
  K2<-max(which(colsbeta!=0))
  
  thres<-c() ## boundaries for abs(beta)
  S0<-list() # estimated relevant sets
  S0c<-list() # estimated irrelevant sets
  
  # the initial values; 
  thres_tmp<-max(abs(BETA[,1]))
  S0_tmp<-which(abs(BETA[,1])>thres_tmp)
  S0c_tmp<-which(abs(BETA[,1])<=thres_tmp)
  
  # loop for the update
  if(K2>=1){
    for (k in 1:K2){
      
      beta_abs<-abs(BETA[,k])
      beta_sort<-sort(beta_abs)
 
      # update the initials
      thres[k]<-ifelse(length(S0c_tmp)>=1,max(beta_abs[S0c_tmp]),0)
      S0_tmp<-which(abs(BETA[,k])>thres[k])
      S0c_tmp<-which(abs(BETA[,k])<=thres[k])
      
      S0[[k]]<-S0_tmp
      S0c[[k]]<-S0c_tmp
      
      if(length(S0c[[k]])>=1){
        
        gap<-diff(c(0,beta_sort))
        
        # the distance between current relevant and irrelevant sets
        gap_10<-ifelse(length(S0c_tmp)==p,0,gap[length(S0c_tmp)+1])
        
        # gap for the current irrelevant set: S0c_tmp
        gap_0<-gap[1:length(S0c_tmp)]
        o1<-which.max(gap_0)
        gap_01<-max(gap_0)
        
        gap_02<-ifelse(o1>1,max(gap[1:(o1-1)]),0)
        
        if(gap_10 <= R * gap_01 & gap_01 >= R*gap_02 ){
          
          thres[k]<-ifelse(o1>1,beta_sort[o1-1],0)

          S0_tmp<-which(abs(BETA[,k])>thres[k])
          S0c_tmp<-which(abs(BETA[,k])<=thres[k])
          
          S0[[k]]<-S0_tmp
          S0c[[k]]<-S0c_tmp
          
        }
      }
      
      
    }
  }
  
  
  index<-rep(1,p)
  
  for(i in 1:p){
    if(all(abs(BETA[i,1:K2])<=thres)) index[i]<-0
  }
  
  
  nz<-which(index==1)
  z<-which(index==0)
  beta_SPSP<-rep(0,p)
  
  if(standardize==FALSE){
    intercept<-0
    if(length(nz)>=1){
      Xc<-X[,nz]
      
      if(length(nz)<=n){
        betac<-lm(Y~Xc-1)$coefficients
      }else{
        betac<-solve(t(Xc)%*%Xc+0.001*diag(length(nz)))%*%t(Xc)%*%Y
      }
      
      beta_SPSP[nz]<-betac  
    }
  }else if(standardize==TRUE){
    if(length(nz)>=1){
      Xc<-X[,nz]
      Xc1<-cbind(1,Xc)
      if(length(nz)<n){
        betac1<-(lm(Y~Xc1-1)$coefficients)
        intercept<-betac1[1]
        betac<-betac1[-1]
      }else{
        betac1<-solve(t(Xc1)%*%Xc1+0.001*diag(length(nz)+1))%*%t(Xc1)%*%Y
        intercept<-betac1[1]
        betac<-betac1[-1]
      }
      
      beta_SPSP[nz]<-betac  
    }else{
      intercept<-lm(Y~1)$coefficients[1]
    }
  }else{
    print("Please specify stardardize==TRUE/FALSE.")
  }
  
  
  list(beta_SPSP=beta_SPSP,S0=S0,thres=thres,
       nonzero=nz,zero=z,R=R,intercept=as.numeric(intercept))
}


