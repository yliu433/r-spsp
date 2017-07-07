

M <- 1

N <- 500 # total simulation
sigma0 <- 1 # noise level


require(MASS)
require(glmnet)
require(ncvreg)

set.seed(100)

source(“1-simulation-model-setups.R”) 
source("SPSP.R") # SPSP function
source("extrafn.R")
source("stabselec.R") # stability selection
###

names <- c(outer(c("L","R","A","S","M"),c("FP","FN","ME"),paste,sep=""))

ReM <- matrix(NA,N,9)
colnames(ReM) <- c("CV","GCV","AIC","BIC","EBIC","Stabs","SPSP1","SPSP5","Ora")
ReM <- as.data.frame(ReM)

for(inames in 1:length(names)){
  assign(names[inames],ReM)
}

ii <- 0  # no intercept
stf1 <- stf2 <- FALSE # no need to standardize

for(j in 1:N){
  
  X <- Xdata[(1+(j-1)*n):(j*n),]
  Y <- Ydata[(1+(j-1)*n):(j*n)]
  
  #**********************
  # P1: Lasso
  
  
  # other methods
  
  #10-fold CV
  cvfit <- cv.glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,standardize = stf1)
  beta_cv0 <- coef(cvfit,s=cvfit$lambda.min)
  beta_cv <- beta_cv0[-1]
  
  LFP[j,"CV"] <- sum(which(beta_cv!=0) %in% zero)
  LFN[j,"CV"] <- sum(which(beta_cv==0) %in% nonzero)
  LME[j,"CV"] <- t(beta_cv-beta)%*%cov(X)%*%(beta_cv-beta)/(sigma^2)

  fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1)
  #plot(fl1,xvar = "lambda")
  K <- dim(fl1$beta)[2]
  LBETA <- as.matrix(fl1$beta[,K:1])

  #GCV
  
  #gc <- apply(LBETA,2,function(t){return(LM_GCV(Y,X,t,intercept = ii))})
  
  gc <- apply(LBETA,2,LM_GCV,y=Y,x=X,intercept = ii)
  #which.min(gc)
  
  beta_gcv <- LBETA[,which.min(gc)]
  
  LFP[j,"GCV"] <- sum(which(beta_gcv!=0) %in% zero)
  LFN[j,"GCV"] <- sum(which(beta_gcv==0) %in% nonzero)
  LME[j,"GCV"] <- t(beta_gcv-beta)%*%cov(X)%*%(beta_gcv-beta)/(sigma^2)

  
  #AIC
  
  #ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
  ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
  #identical(ac1,ac)
  #which.min(ac)
  
  beta_aic <- LBETA[,which.min(ac)]
  
  LFP[j,"AIC"] <- sum(which(beta_aic!=0) %in% zero)
  LFN[j,"AIC"] <- sum(which(beta_aic==0) %in% nonzero)
  LME[j,"AIC"] <- t(beta_aic-beta)%*%cov(X)%*%(beta_aic-beta)/(sigma^2)

  #BIC
  
  #bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
  bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
  #identical(bc1,bc)
  #which.min(bc)
  
  beta_bic <- LBETA[,which.min(bc)]
  
  LFP[j,"BIC"] <- sum(which(beta_bic!=0) %in% zero)
  LFN[j,"BIC"] <- sum(which(beta_bic==0) %in% nonzero)
  LME[j,"BIC"] <- t(beta_bic-beta)%*%cov(X)%*%(beta_bic-beta)/(sigma^2)

  #EBIC
  
  #ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
  ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
  #identical(ebc1,ebc)
  #which.min(ebc)
  
  beta_ebic <- LBETA[,which.min(ebc)]
  
  LFP[j,"EBIC"] <- sum(which(beta_ebic!=0) %in% zero)
  LFN[j,"EBIC"] <- sum(which(beta_ebic==0) %in% nonzero)
  LME[j,"EBIC"] <- t(beta_ebic-beta)%*%cov(X)%*%(beta_ebic-beta)/(sigma^2)


  
  
  #Oracle
  
  FP_tmp <- apply(LBETA,2,function(x){sum(x[zero]!=0)})
  FN_tmp <- apply(LBETA,2,function(x){sum(x[nonzero]==0)})
  PS_tmp <- FP_tmp+FN_tmp
  #o_tmp <- max(which(PS_tmp==min(PS_tmp)))
  
  beta_ora <- LBETA[,which.min(PS_tmp)]
  
  LFP[j,"Ora"] <- sum(which(beta_ora!=0) %in% zero)
  LFN[j,"Ora"] <- sum(which(beta_ora==0) %in% nonzero)
  LME[j,"Ora"] <- t(beta_ora-beta)%*%cov(X)%*%(beta_ora-beta)/(sigma^2)
  
  
  # stab selection
  
  stabss <- stab_lasso(X,Y)
  LFP[j,"Stabs"] <- sum(stabss$nonzero %in% zero)
  LFN[j,"Stabs"] <- sum(stabss$zero %in% nonzero)
  beta_stabs <- stabss$beta_stab
  LME[j,"Stabs"] <- t(beta_stabs-beta)%*%cov(X)%*%(beta_stabs-beta)/(sigma^2)
  
  # SPSP solution path
  fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
                standardize = stf1,lambda.min.ratio = 0.01)
  #plot(fl1,xvar = "lambda")
  K <- dim(fl1$beta)[2]
  LBETA <- as.matrix(fl1$beta[,K:1])
  
  # SPSP
  r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
  LFP[j,"SPSP1"] <- sum(r_spsp1$nonzero %in% zero)
  LFN[j,"SPSP1"] <- sum(r_spsp1$zero %in% nonzero)
  beta_spsp1 <- r_spsp1$beta_SPSP
  LME[j,"SPSP1"] <- t(beta_spsp1-beta)%*%cov(X)%*%(beta_spsp1-beta)/(sigma^2)
  
  r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
  LFP[j,"SPSP5"] <- sum(r_spsp5$nonzero %in% zero)
  LFN[j,"SPSP5"] <- sum(r_spsp5$zero %in% nonzero)
  beta_spsp5 <- r_spsp5$beta_SPSP
  LME[j,"SPSP5"] <- t(beta_spsp5-beta)%*%cov(X)%*%(beta_spsp5-beta)/(sigma^2)
  
  
  #**********************
  # P2: Ridge
  
  # using the same grid as lasso
  fl0 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,lambda.min.ratio = 0.01)
  
  fl1 <- glmnet(X,Y,family="gaussian",alpha=0,intercept=ii,standardize = stf1,lambda=fl0$lambda)
  #plot(fr)
  #print(fr)
  #plot(fl1,xvar = "lambda")
  K <- dim(fl1$beta)[2]
  LBETA <- as.matrix(fl1$beta[,K:1])
  
  # SPSP
  r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
  RFP[j,"SPSP1"] <- sum(r_spsp1$nonzero %in% zero)
  RFN[j,"SPSP1"] <- sum(r_spsp1$zero %in% nonzero)
  beta_spsp1 <- r_spsp1$beta_SPSP
  RME[j,"SPSP1"] <- t(beta_spsp1-beta)%*%cov(X)%*%(beta_spsp1-beta)/(sigma^2)
  
  r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
  RFP[j,"SPSP5"] <- sum(r_spsp5$nonzero %in% zero)
  RFN[j,"SPSP5"] <- sum(r_spsp5$zero %in% nonzero)
  beta_spsp5 <- r_spsp5$beta_SPSP
  RME[j,"SPSP5"] <- t(beta_spsp5-beta)%*%cov(X)%*%(beta_spsp5-beta)/(sigma^2)
  
  
  #**********************
  # P3: adaLasso
  
  
  # solution path
  cv_ridge  <-  cv.glmnet(X, Y, family='gaussian', alpha=0, standardize=stf1,intercept=ii)
  w3  <-  1/abs(matrix(coef(cv_ridge, s=cv_ridge$lambda.min)[-1]))^1
  
  
  fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,penalty.factor = w3,lambda.min.ratio = 0.01)
  #plot(fl1,xvar = "lambda")
  K <- dim(fl1$beta)[2]
  LBETA <- as.matrix(fl1$beta[,K:1])
  
  # SPSP
  r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
  AFP[j,"SPSP1"] <- sum(r_spsp1$nonzero %in% zero)
  AFN[j,"SPSP1"] <- sum(r_spsp1$zero %in% nonzero)
  beta_spsp1 <- r_spsp1$beta_SPSP
  AME[j,"SPSP1"] <- t(beta_spsp1-beta)%*%cov(X)%*%(beta_spsp1-beta)/(sigma^2)
  
  r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
  AFP[j,"SPSP5"] <- sum(r_spsp5$nonzero %in% zero)
  AFN[j,"SPSP5"] <- sum(r_spsp5$zero %in% nonzero)
  beta_spsp5 <- r_spsp5$beta_SPSP
  AME[j,"SPSP5"] <- t(beta_spsp5-beta)%*%cov(X)%*%(beta_spsp5-beta)/(sigma^2)
  
  # other methods
  
  # stab selection
  
  stabss <- stab_adalasso(X,Y)
  AFP[j,"Stabs"] <- sum(stabss$nonzero %in% zero)
  AFN[j,"Stabs"] <- sum(stabss$zero %in% nonzero)
  beta_stabs <- stabss$beta_stab
  AME[j,"Stabs"] <- t(beta_stabs-beta)%*%cov(X)%*%(beta_stabs-beta)/(sigma^2)
  
  #10-fold CV
  cvfit<-cv.glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,standardize = stf1,
                   penalty.factor = w3)
  beta_cv0<-coef(cvfit,s=cvfit$lambda.min)
  beta_cv<-beta_cv0[-1]
  
  AFP[j,"CV"] <- sum(which(beta_cv!=0) %in% zero)
  AFN[j,"CV"] <- sum(which(beta_cv==0) %in% nonzero)
  AME[j,"CV"] <- t(beta_cv-beta)%*%cov(X)%*%(beta_cv-beta)/(sigma^2)
  
  
  #GCV
  
  #gc <- apply(LBETA,2,function(t){return(LM_GCV(Y,X,t,intercept = ii))})
  fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,penalty.factor = w3)
  #plot(fl1,xvar = "lambda")
  K <- dim(fl1$beta)[2]
  LBETA <- as.matrix(fl1$beta[,K:1])

  gc <- apply(LBETA,2,LM_GCV,y=Y,x=X,intercept = ii)
  #which.min(gc)
  
  beta_gcv <- LBETA[,which.min(gc)]
  
  AFP[j,"GCV"] <- sum(which(beta_gcv!=0) %in% zero)
  AFN[j,"GCV"] <- sum(which(beta_gcv==0) %in% nonzero)
  AME[j,"GCV"] <- t(beta_gcv-beta)%*%cov(X)%*%(beta_gcv-beta)/(sigma^2)
  
  
  #AIC
  
  #ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
  ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
  #identical(ac1,ac)
  #which.min(ac)
  
  beta_aic <- LBETA[,which.min(ac)]
  
  AFP[j,"AIC"] <- sum(which(beta_aic!=0) %in% zero)
  AFN[j,"AIC"] <- sum(which(beta_aic==0) %in% nonzero)
  AME[j,"AIC"] <- t(beta_aic-beta)%*%cov(X)%*%(beta_aic-beta)/(sigma^2)
  
  #BIC
  
  #bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
  bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
  #identical(bc1,bc)
  #which.min(bc)
  
  beta_bic <- LBETA[,which.min(bc)]
  
  AFP[j,"BIC"] <- sum(which(beta_bic!=0) %in% zero)
  AFN[j,"BIC"] <- sum(which(beta_bic==0) %in% nonzero)
  AME[j,"BIC"] <- t(beta_bic-beta)%*%cov(X)%*%(beta_bic-beta)/(sigma^2)
  
  #EBIC
  
  #ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
  ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
  #identical(ebc1,ebc)
  #which.min(ebc)
  
  beta_ebic <- LBETA[,which.min(ebc)]
  
  AFP[j,"EBIC"] <- sum(which(beta_ebic!=0) %in% zero)
  AFN[j,"EBIC"] <- sum(which(beta_ebic==0) %in% nonzero)
  AME[j,"EBIC"] <- t(beta_ebic-beta)%*%cov(X)%*%(beta_ebic-beta)/(sigma^2)
  
  
  #Oracle
  
  FP_tmp <- apply(LBETA,2,function(x){sum(x[zero]!=0)})
  FN_tmp <- apply(LBETA,2,function(x){sum(x[nonzero]==0)})
  PS_tmp <- FP_tmp+FN_tmp
  #o_tmp <- max(which(PS_tmp==min(PS_tmp)))
  
  beta_ora <- LBETA[,which.min(PS_tmp)]
  
  AFP[j,"Ora"] <- sum(which(beta_ora!=0) %in% zero)
  AFN[j,"Ora"] <- sum(which(beta_ora==0) %in% nonzero)
  AME[j,"Ora"] <- t(beta_ora-beta)%*%cov(X)%*%(beta_ora-beta)/(sigma^2)
  
  
  #********************
  #  P4: SCAD
  
  fs <- ncvreg(X,Y,family = "gaussian",penalty="SCAD",lambda.min = 0.1)
  #plot(fs)
  #plot(fl1,xvar = "lambda")
  K <- dim(fs$beta)[2]
  LBETA <- as.matrix(fs$beta[-1,K:1])
  
  # SPSP
  r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
  SFP[j,"SPSP1"] <- sum(r_spsp1$nonzero %in% zero)
  SFN[j,"SPSP1"] <- sum(r_spsp1$zero %in% nonzero)
  beta_spsp1 <- r_spsp1$beta_SPSP
  SME[j,"SPSP1"] <- t(beta_spsp1-beta)%*%cov(X)%*%(beta_spsp1-beta)/(sigma^2)
  
  r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
  SFP[j,"SPSP5"] <- sum(r_spsp5$nonzero %in% zero)
  SFN[j,"SPSP5"] <- sum(r_spsp5$zero %in% nonzero)
  beta_spsp5 <- r_spsp5$beta_SPSP
  SME[j,"SPSP5"] <- t(beta_spsp5-beta)%*%cov(X)%*%(beta_spsp5-beta)/(sigma^2)
  
  # other methods
  # stab selection
  
  stabss <- stab_scad(X,Y)
  SFP[j,"Stabs"] <- sum(stabss$nonzero %in% zero)
  SFN[j,"Stabs"] <- sum(stabss$zero %in% nonzero)
  beta_stabs <- stabss$beta_stab
  SME[j,"Stabs"] <- t(beta_stabs-beta)%*%cov(X)%*%(beta_stabs-beta)/(sigma^2)
  #10-fold CV
  cvfit<-cv.ncvreg(X,Y,family = "gaussian",penalty="SCAD")
  beta_cv0<-fs$beta[,cvfit$min]
  beta_cv<-beta_cv0[-1]
  
  SFP[j,"CV"] <- sum(which(beta_cv!=0) %in% zero)
  SFN[j,"CV"] <- sum(which(beta_cv==0) %in% nonzero)
  SME[j,"CV"] <- t(beta_cv-beta)%*%cov(X)%*%(beta_cv-beta)/(sigma^2)
  
  
  #GCV
  fs <- ncvreg(X,Y,family = "gaussian",penalty="SCAD")
  #plot(fs)
  #plot(fl1,xvar = "lambda")
  K <- dim(fs$beta)[2]
  LBETA <- as.matrix(fs$beta[-1,K:1])

  #gc <- apply(LBETA,2,function(t){return(LM_GCV(Y,X,t,intercept = ii))})
  
  gc <- apply(LBETA,2,LM_GCV,y=Y,x=X,intercept = ii)
  #which.min(gc)
  
  beta_gcv <- LBETA[,which.min(gc)]
  
  SFP[j,"GCV"] <- sum(which(beta_gcv!=0) %in% zero)
  SFN[j,"GCV"] <- sum(which(beta_gcv==0) %in% nonzero)
  SME[j,"GCV"] <- t(beta_gcv-beta)%*%cov(X)%*%(beta_gcv-beta)/(sigma^2)
  
  
  #AIC
  
  #ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
  ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
  #identical(ac1,ac)
  #which.min(ac)
  
  beta_aic <- LBETA[,which.min(ac)]
  
  SFP[j,"AIC"] <- sum(which(beta_aic!=0) %in% zero)
  SFN[j,"AIC"] <- sum(which(beta_aic==0) %in% nonzero)
  SME[j,"AIC"] <- t(beta_aic-beta)%*%cov(X)%*%(beta_aic-beta)/(sigma^2)
  
  #BIC
  
  #bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
  bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
  #identical(bc1,bc)
  #which.min(bc)
  
  beta_bic <- LBETA[,which.min(bc)]
  
  SFP[j,"BIC"] <- sum(which(beta_bic!=0) %in% zero)
  SFN[j,"BIC"] <- sum(which(beta_bic==0) %in% nonzero)
  SME[j,"BIC"] <- t(beta_bic-beta)%*%cov(X)%*%(beta_bic-beta)/(sigma^2)
  
  #EBIC
  
  #ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
  ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
  #identical(ebc1,ebc)
  #which.min(ebc)
  
  beta_ebic <- LBETA[,which.min(ebc)]
  
  SFP[j,"EBIC"] <- sum(which(beta_ebic!=0) %in% zero)
  SFN[j,"EBIC"] <- sum(which(beta_ebic==0) %in% nonzero)
  SME[j,"EBIC"] <- t(beta_ebic-beta)%*%cov(X)%*%(beta_ebic-beta)/(sigma^2)
  
  
  #Oracle
  
  FP_tmp <- apply(LBETA,2,function(x){sum(x[zero]!=0)})
  FN_tmp <- apply(LBETA,2,function(x){sum(x[nonzero]==0)})
  PS_tmp <- FP_tmp+FN_tmp
  #o_tmp <- max(which(PS_tmp==min(PS_tmp)))
  
  beta_ora <- LBETA[,which.min(PS_tmp)]
  
  SFP[j,"Ora"] <- sum(which(beta_ora!=0) %in% zero)
  SFN[j,"Ora"] <- sum(which(beta_ora==0) %in% nonzero)
  SME[j,"Ora"] <- t(beta_ora-beta)%*%cov(X)%*%(beta_ora-beta)/(sigma^2)
  
  #********************
  #  P5: MCP
  
  fs <- ncvreg(X,Y,family = "gaussian",penalty="MCP",lambda.min = 0.1)
  #plot(fs)
  K <- dim(fs$beta)[2]
  LBETA <- as.matrix(fs$beta[-1,K:1])
  
  # SPSP
  r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
  MFP[j,"SPSP1"] <- sum(r_spsp1$nonzero %in% zero)
  MFN[j,"SPSP1"] <- sum(r_spsp1$zero %in% nonzero)
  beta_spsp1 <- r_spsp1$beta_SPSP
  MME[j,"SPSP1"] <- t(beta_spsp1-beta)%*%cov(X)%*%(beta_spsp1-beta)/(sigma^2)
  
  r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
  MFP[j,"SPSP5"] <- sum(r_spsp5$nonzero %in% zero)
  MFN[j,"SPSP5"] <- sum(r_spsp5$zero %in% nonzero)
  beta_spsp5 <- r_spsp5$beta_SPSP
  MME[j,"SPSP5"] <- t(beta_spsp5-beta)%*%cov(X)%*%(beta_spsp5-beta)/(sigma^2)
  
  # other methods
  # stab selection
  
  stabss <- stab_mcp(X,Y)
  MFP[j,"Stabs"] <- sum(stabss$nonzero %in% zero)
  MFN[j,"Stabs"] <- sum(stabss$zero %in% nonzero)
  beta_stabs <- stabss$beta_stab
  MME[j,"Stabs"] <- t(beta_stabs-beta)%*%cov(X)%*%(beta_stabs-beta)/(sigma^2)
  
  #10-fold CV
  cvfit<-cv.ncvreg(X,Y,family = "gaussian",penalty="MCP")
  beta_cv0<-fs$beta[,cvfit$min]
  beta_cv<-beta_cv0[-1]
  
  MFP[j,"CV"] <- sum(which(beta_cv!=0) %in% zero)
  MFN[j,"CV"] <- sum(which(beta_cv==0) %in% nonzero)
  MME[j,"CV"] <- t(beta_cv-beta)%*%cov(X)%*%(beta_cv-beta)/(sigma^2)
  
  
  #GCV
  fs <- ncvreg(X,Y,family = "gaussian",penalty="MCP")
  #plot(fs)
  K <- dim(fs$beta)[2]
  LBETA <- as.matrix(fs$beta[-1,K:1])

  #gc <- apply(LBETA,2,function(t){return(LM_GCV(Y,X,t,intercept = ii))})
  
  gc <- apply(LBETA,2,LM_GCV,y=Y,x=X,intercept = ii)
  #which.min(gc)
  
  beta_gcv <- LBETA[,which.min(gc)]
  
  MFP[j,"GCV"] <- sum(which(beta_gcv!=0) %in% zero)
  MFN[j,"GCV"] <- sum(which(beta_gcv==0) %in% nonzero)
  MME[j,"GCV"] <- t(beta_gcv-beta)%*%cov(X)%*%(beta_gcv-beta)/(sigma^2)
  
  
  #AIC
  
  #ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
  ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
  #identical(ac1,ac)
  #which.min(ac)
  
  beta_aic <- LBETA[,which.min(ac)]
  
  MFP[j,"AIC"] <- sum(which(beta_aic!=0) %in% zero)
  MFN[j,"AIC"] <- sum(which(beta_aic==0) %in% nonzero)
  MME[j,"AIC"] <- t(beta_aic-beta)%*%cov(X)%*%(beta_aic-beta)/(sigma^2)
  
  #BIC
  
  #bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
  bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
  #identical(bc1,bc)
  #which.min(bc)
  
  beta_bic <- LBETA[,which.min(bc)]
  
  MFP[j,"BIC"] <- sum(which(beta_bic!=0) %in% zero)
  MFN[j,"BIC"] <- sum(which(beta_bic==0) %in% nonzero)
  MME[j,"BIC"] <- t(beta_bic-beta)%*%cov(X)%*%(beta_bic-beta)/(sigma^2)
  
  #EBIC
  
  #ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
  ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
  #identical(ebc1,ebc)
  #which.min(ebc)
  
  beta_ebic <- LBETA[,which.min(ebc)]
  
  MFP[j,"EBIC"] <- sum(which(beta_ebic!=0) %in% zero)
  MFN[j,"EBIC"] <- sum(which(beta_ebic==0) %in% nonzero)
  MME[j,"EBIC"] <- t(beta_ebic-beta)%*%cov(X)%*%(beta_ebic-beta)/(sigma^2)
  
  
  #Oracle
  
  FP_tmp <- apply(LBETA,2,function(x){sum(x[zero]!=0)})
  FN_tmp <- apply(LBETA,2,function(x){sum(x[nonzero]==0)})
  PS_tmp <- FP_tmp+FN_tmp
  #o_tmp <- max(which(PS_tmp==min(PS_tmp)))
  
  beta_ora <- LBETA[,which.min(PS_tmp)]
  
  MFP[j,"Ora"] <- sum(which(beta_ora!=0) %in% zero)
  MFN[j,"Ora"] <- sum(which(beta_ora==0) %in% nonzero)
  MME[j,"Ora"] <- t(beta_ora-beta)%*%cov(X)%*%(beta_ora-beta)/(sigma^2)
  
  
  if(j%%(N/5)==1){
    cat(j,"=>")
  }
}

#######################
# mean summary

# lasso
rp <- round(apply(LFP,2,mean,na.rm=T),3);rp2 <- round(apply(LFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(LFN,2,mean,na.rm=T),3);rn2 <- round(apply(LFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(LME,2,mean,na.rm=T),3);re2 <- round(apply(LME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f1 <- rbind(rp,rp2,rn,rn2,re,re2)

# alasso
rp <- round(apply(AFP,2,mean,na.rm=T),3);rp2 <- round(apply(AFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(AFN,2,mean,na.rm=T),3);rn2 <- round(apply(AFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(AME,2,mean,na.rm=T),3);re2 <- round(apply(AME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f2 <- rbind(rp,rp2,rn,rn2,re,re2)

# scad
rp <- round(apply(SFP,2,mean,na.rm=T),3);rp2 <- round(apply(SFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(SFN,2,mean,na.rm=T),3);rn2 <- round(apply(SFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(SME,2,mean,na.rm=T),3);re2 <- round(apply(SME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f3 <- rbind(rp,rp2,rn,rn2,re,re2)

# mcp
rp <- round(apply(MFP,2,mean,na.rm=T),3);rp2 <- round(apply(MFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(MFN,2,mean,na.rm=T),3);rn2 <- round(apply(MFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(MME,2,mean,na.rm=T),3);re2 <- round(apply(MME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f4 <- rbind(rp,rp2,rn,rn2,re,re2)

# ridge
rp <- round(apply(RFP,2,mean,na.rm=T),3);rp2 <- round(apply(RFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(RFN,2,mean,na.rm=T),3);rn2 <- round(apply(RFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(RME,2,mean,na.rm=T),3);re2 <- round(apply(RME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f5 <- rbind(rp,rp2,rn,rn2,re,re2)[,c("SPSP1","SPSP5")]

################################

FR <- rbind(f1,f2,f3,f4)
  
c1<-c("Lasso",rep("",5),"adaLasso",rep("",5),"SCAD",rep("",5),"MCP",rep("",5))
c2<-rep(c("FP","","FN","","ME",""),4)
FR<-cbind(c1,c2,FR)
rownames(FR)<-NULL

ind <- apply(FR,1,function(x) !any(grepl("\\(",x)))


cat("\n--------------------\n")
cat("The mean summary:\n")
print(FR[ind,],quote=FALSE)
cat("Ridge:\n")
print(f5[c(1,3,5),],quote=FALSE)

FR1 <- FR;f51 <- f5

#######################
# median summary

# lasso
rp <- round(apply(LFP,2,median,na.rm=T),3);rp2 <- round(1.253*apply(LFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(LFN,2,median,na.rm=T),3);rn2 <- round(1.253*apply(LFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(LME,2,median,na.rm=T),3);re2 <- round(1.253*apply(LME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f1 <- rbind(rp,rp2,rn,rn2,re,re2)

# alasso
rp <- round(apply(AFP,2,median,na.rm=T),3);rp2 <- round(1.253*apply(AFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(AFN,2,median,na.rm=T),3);rn2 <- round(1.253*apply(AFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(AME,2,median,na.rm=T),3);re2 <- round(1.253*apply(AME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f2 <- rbind(rp,rp2,rn,rn2,re,re2)

# scad
rp <- round(apply(SFP,2,median,na.rm=T),3);rp2 <- round(1.253*apply(SFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(SFN,2,median,na.rm=T),3);rn2 <- round(1.253*apply(SFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(SME,2,median,na.rm=T),3);re2 <- round(1.253*apply(SME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f3 <- rbind(rp,rp2,rn,rn2,re,re2)

# mcp
rp <- round(apply(MFP,2,median,na.rm=T),3);rp2 <- round(1.253*apply(MFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(MFN,2,median,na.rm=T),3);rn2 <- round(1.253*apply(MFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(MME,2,median,na.rm=T),3);re2 <- round(1.253*apply(MME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f4 <- rbind(rp,rp2,rn,rn2,re,re2)

# ridge
rp <- round(apply(RFP,2,median,na.rm=T),3);rp2 <- round(1.253*apply(RFP,2,sd,na.rm=T)/sqrt(N),3)
rn <- round(apply(RFN,2,median,na.rm=T),3);rn2 <- round(1.253*apply(RFN,2,sd,na.rm=T)/sqrt(N),3)
re <- round(apply(RME,2,median,na.rm=T),3);re2 <- round(1.253*apply(RME,2,sd,na.rm=T)/sqrt(N),3)

rp2 <- paste0("(", rp2,")")
rn2 <- paste0("(", rn2,")")
re2 <- paste0("(", re2,")")

f5 <- rbind(rp,rp2,rn,rn2,re,re2)[,c("SPSP1","SPSP5")]

################################

FR <- rbind(f1,f2,f3,f4)

c1<-c("Lasso",rep("",5),"adaLasso",rep("",5),"SCAD",rep("",5),"MCP",rep("",5))
c2<-rep(c("FP","","FN","","ME",""),4)
FR<-cbind(c1,c2,FR)
rownames(FR)<-NULL

ind <- apply(FR,1,function(x) !any(grepl("\\(",x)))


cat("\n--------------------\n")
cat("The median summary:\n")
print(FR[ind,],quote=FALSE)
cat("Ridge:\n")
print(f5[c(1,3,5),],quote=FALSE)

FR2 <- FR;f52 <- f5

q(save="no")

