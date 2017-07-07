
require(MASS)
require(glmnet)
require(ncvreg)
source("SPSP.R") # SPSP function
source("extrafn.R")
source("stabselec.R") # stability selection

load("gbm.Rdata")

X <- Xs1000
Y <- Yc

RM <- matrix(NA,10,8)
colnames(RM) <- c("CV","GCV","AIC","BIC","EBIC","Stabs","SPSP1","SPSP5")
rownames(RM) <- paste(rep(c("L","A","S","M","R"),each=2),rep(c("Num","PME"),5),sep="")


ii <- 0  # no intercept
stf1 <- stf2 <- FALSE # no need to standardize

set.seed(100)

#**********************
# P1: Lasso


# other methods

#10-fold CV
cvfit <- cv.glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,standardize = stf1)
beta_cv0 <- coef(cvfit,s=cvfit$lambda.min)
beta_cv <- beta_cv0[-1]

RM["LNum","CV"] <- sum(beta_cv!=0)

RM["LPME","CV"] <- mean((Y_test-X_test%*%beta_cv)^2)

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

RM["LNum","GCV"] <- sum(beta_gcv!=0)

RM["LPME","GCV"] <- mean((Y_test-X_test%*%beta_gcv)^2)


#AIC

#ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
#identical(ac1,ac)
#which.min(ac)

beta_aic <- LBETA[,which.min(ac)]

RM["LNum","AIC"] <- sum(beta_aic!=0)

RM["LPME","AIC"] <- mean((Y_test-X_test%*%beta_aic)^2)

#BIC

#bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
#identical(bc1,bc)
#which.min(bc)

beta_bic <- LBETA[,which.min(bc)]

RM["LNum","BIC"] <- sum(beta_bic!=0)

RM["LPME","BIC"] <- mean((Y_test-X_test%*%beta_bic)^2)

#EBIC

#ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
#identical(ebc1,ebc)
#which.min(ebc)

beta_ebic <- LBETA[,which.min(ebc)]

RM["LNum","EBIC"] <- sum(beta_ebic!=0)

RM["LPME","EBIC"] <- mean((Y_test-X_test%*%beta_ebic)^2)





# stab selection

stabss <- stab_lasso(X,Y)
beta_stabs <- stabss$beta_stab

RM["LNum","Stabs"] <- sum(beta_stabs!=0)

RM["LPME","Stabs"] <- mean((Y_test-X_test%*%beta_stabs)^2)

# SPSP solution path
fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,lambda.min.ratio = 0.05)
#plot(fl1,xvar = "lambda")
K <- dim(fl1$beta)[2]
LBETA <- as.matrix(fl1$beta[,K:1])

# SPSP
r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
beta_spsp1 <- r_spsp1$beta_SPSP

RM["LNum","SPSP1"] <- sum(beta_spsp1!=0)

RM["LPME","SPSP1"] <- mean((Y_test-X_test%*%beta_spsp1)^2)


r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
beta_spsp5 <- r_spsp5$beta_SPSP

RM["LNum","SPSP5"] <- sum(beta_spsp5!=0)

RM["LPME","SPSP5"] <- mean((Y_test-X_test%*%beta_spsp5)^2)


#**********************
# P2: Ridge

# using the same grid as lasso
fl0 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,lambda.min.ratio = 0.05)

fl1 <- glmnet(X,Y,family="gaussian",alpha=0,intercept=ii,standardize = stf1,lambda=fl0$lambda)
#plot(fr)
#print(fr)
#plot(fl1,xvar = "lambda")
K <- dim(fl1$beta)[2]
LBETA <- as.matrix(fl1$beta[,K:1])

r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
beta_spsp1 <- r_spsp1$beta_SPSP

RM["RNum","SPSP1"] <- sum(beta_spsp1!=0)

RM["RPME","SPSP1"] <- mean((Y_test-X_test%*%beta_spsp1)^2)


r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
beta_spsp5 <- r_spsp5$beta_SPSP

RM["RNum","SPSP5"] <- sum(beta_spsp5!=0)

RM["RPME","SPSP5"] <- mean((Y_test-X_test%*%beta_spsp5)^2)

#**********************
# P3: adaLasso


# solution path
cv_ridge  <-  cv.glmnet(X, Y, family='gaussian', alpha=0, standardize=stf1,intercept=ii)
w3  <-  1/abs(matrix(coef(cv_ridge, s=cv_ridge$lambda.min)[-1]))^1


fl1 <- glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,
              standardize = stf1,penalty.factor = w3,lambda.min.ratio = 0.05)
#plot(fl1,xvar = "lambda")
K <- dim(fl1$beta)[2]
LBETA <- as.matrix(fl1$beta[,K:1])

# SPSP
r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
beta_spsp1 <- r_spsp1$beta_SPSP
RM["ANum","SPSP1"] <- sum(beta_spsp1!=0)

RM["APME","SPSP1"] <- mean((Y_test-X_test%*%beta_spsp1)^2)

r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
beta_spsp5 <- r_spsp5$beta_SPSP
RM["ANum","SPSP5"] <- sum(beta_spsp5!=0)

RM["APME","SPSP5"] <- mean((Y_test-X_test%*%beta_spsp5)^2)

# other methods

# stab selection

stabss <- stab_adalasso(X,Y)
beta_stabs <- stabss$beta_stab

RM["ANum","Stabs"] <- sum(beta_stabs!=0)

RM["APME","Stabs"] <- mean((Y_test-X_test%*%beta_stabs)^2)

#10-fold CV
cvfit<-cv.glmnet(X,Y,family="gaussian",alpha=1,intercept=ii,standardize = stf1,
                 penalty.factor = w3)
beta_cv0<-coef(cvfit,s=cvfit$lambda.min)
beta_cv<-beta_cv0[-1]

RM["ANum","CV"] <- sum(beta_cv!=0)

RM["APME","CV"] <- mean((Y_test-X_test%*%beta_cv)^2)


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

RM["ANum","GCV"] <- sum(beta_gcv!=0)

RM["APME","GCV"] <- mean((Y_test-X_test%*%beta_gcv)^2)


#AIC

#ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
#identical(ac1,ac)
#which.min(ac)

beta_aic <- LBETA[,which.min(ac)]

RM["ANum","AIC"] <- sum(beta_aic!=0)

RM["APME","AIC"] <- mean((Y_test-X_test%*%beta_aic)^2)

#BIC

#bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
#identical(bc1,bc)
#which.min(bc)

beta_bic <- LBETA[,which.min(bc)]

RM["ANum","BIC"] <- sum(beta_bic!=0)

RM["APME","BIC"] <- mean((Y_test-X_test%*%beta_bic)^2)

#EBIC

#ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
#identical(ebc1,ebc)
#which.min(ebc)

beta_ebic <- LBETA[,which.min(ebc)]

RM["ANum","EBIC"] <- sum(beta_ebic!=0)

RM["APME","EBIC"] <- mean((Y_test-X_test%*%beta_ebic)^2)

#********************
#  P4: SCAD

fs <- ncvreg(X,Y,family = "gaussian",penalty="SCAD",lambda.min = 0.1)
#plot(fs)
#plot(fl1,xvar = "lambda")
K <- dim(fs$beta)[2]
LBETA <- as.matrix(fs$beta[-1,K:1])

# SPSP
r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
beta_spsp1 <- r_spsp1$beta_SPSP

RM["SNum","SPSP1"] <- sum(beta_spsp1!=0)

RM["SPME","SPSP1"] <- mean((Y_test-X_test%*%beta_spsp1)^2)

r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
beta_spsp5 <- r_spsp5$beta_SPSP

RM["SNum","SPSP5"] <- sum(beta_spsp5!=0)

RM["SPME","SPSP5"] <- mean((Y_test-X_test%*%beta_spsp5)^2)
# other methods
# stab selection

stabss <- stab_scad(X,Y)
beta_stabs <- stabss$beta_stab
RM["SNum","Stabs"] <- sum(beta_stabs!=0)

RM["SPME","Stabs"] <- mean((Y_test-X_test%*%beta_stabs)^2)

#10-fold CV
cvfit<-cv.ncvreg(X,Y,family = "gaussian",penalty="SCAD")
beta_cv0<-fs$beta[,cvfit$min]
beta_cv<-beta_cv0[-1]

RM["SNum","CV"] <- sum(beta_cv!=0)

RM["SPME","CV"] <- mean((Y_test-X_test%*%beta_cv)^2)


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

RM["SNum","GCV"] <- sum(beta_gcv!=0)

RM["SPME","GCV"] <- mean((Y_test-X_test%*%beta_gcv)^2)


#AIC

#ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
#identical(ac1,ac)
#which.min(ac)

beta_aic <- LBETA[,which.min(ac)]

RM["SNum","AIC"] <- sum(beta_aic!=0)

RM["SPME","AIC"] <- mean((Y_test-X_test%*%beta_aic)^2)

#BIC

#bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
#identical(bc1,bc)
#which.min(bc)

beta_bic <- LBETA[,which.min(bc)]

RM["SNum","BIC"] <- sum(beta_bic!=0)

RM["SPME","BIC"] <- mean((Y_test-X_test%*%beta_bic)^2)

#EBIC

#ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
#identical(ebc1,ebc)
#which.min(ebc)

beta_ebic <- LBETA[,which.min(ebc)]

RM["SNum","EBIC"] <- sum(beta_ebic!=0)

RM["SPME","EBIC"] <- mean((Y_test-X_test%*%beta_ebic)^2)



#********************
#  P5: MCP

fs <- ncvreg(X,Y,family = "gaussian",penalty="MCP",lambda.min = 0.1)
#plot(fs)
K <- dim(fs$beta)[2]
LBETA <- as.matrix(fs$beta[-1,K:1])

# SPSP
r_spsp1 <- SPSP(Y,X,LBETA,init = 1,standardize = stf2)
beta_spsp1 <- r_spsp1$beta_SPSP

RM["MNum","SPSP1"] <- sum(beta_spsp1!=0)

RM["MPME","SPSP1"] <- mean((Y_test-X_test%*%beta_spsp1)^2)

r_spsp5 <- SPSP(Y,X,LBETA,init = 5,standardize = stf2)
beta_spsp5 <- r_spsp5$beta_SPSP

RM["MNum","SPSP5"] <- sum(beta_spsp5!=0)

RM["MPME","SPSP5"] <- mean((Y_test-X_test%*%beta_spsp5)^2)

# other methods
# stab selection

stabss <- stab_mcp(X,Y)
beta_stabs <- stabss$beta_stab
RM["MNum","Stabs"] <- sum(beta_stabs!=0)

RM["MPME","Stabs"] <- mean((Y_test-X_test%*%beta_stabs)^2)

#10-fold CV
cvfit<-cv.ncvreg(X,Y,family = "gaussian",penalty="MCP")
beta_cv0<-fs$beta[,cvfit$min]
beta_cv<-beta_cv0[-1]

RM["MNum","CV"] <- sum(beta_cv!=0)

RM["MPME","CV"] <- mean((Y_test-X_test%*%beta_cv)^2)


#GCV
fs <- ncvreg(X,Y,family = "gaussian",penalty="MCP")
#plot(fs)
K <- dim(fs$beta)[2]
LBETA <- as.matrix(fs$beta[-1,K:1])

#gc <- apply(LBETA,2,function(t){return(LM_GCV(Y,X,t,intercept = ii))})

gc <- apply(LBETA,2,LM_GCV,y=Y,x=X,intercept = ii)
#which.min(gc)

beta_gcv <- LBETA[,which.min(gc)]

RM["MNum","GCV"] <- sum(beta_gcv!=0)

RM["MPME","GCV"] <- mean((Y_test-X_test%*%beta_gcv)^2)


#AIC

#ac1 <- apply(LBETA,2,function(t){return(LM_AIC(Y,X,t,intercept = ii))})
ac <- apply(LBETA,2,LM_AIC,y=Y,x=X,intercept = ii)
#identical(ac1,ac)
#which.min(ac)

beta_aic <- LBETA[,which.min(ac)]

RM["MNum","AIC"] <- sum(beta_aic!=0)

RM["MPME","AIC"] <- mean((Y_test-X_test%*%beta_aic)^2)

#BIC

#bc1 <- apply(LBETA,2,function(t){return(LM_BIC(Y,X,t,intercept = ii))})
bc <- apply(LBETA,2,LM_BIC,y=Y,x=X,intercept = ii)
#identical(bc1,bc)
#which.min(bc)

beta_bic <- LBETA[,which.min(bc)]

RM["MNum","BIC"] <- sum(beta_bic!=0)

RM["MPME","BIC"] <- mean((Y_test-X_test%*%beta_bic)^2)

#EBIC

#ebc1 <- apply(LBETA,2,function(t){return(LM_EBIC(Y,X,t,intercept = ii,gamma=1))})
ebc <- apply(LBETA,2,LM_EBIC,y=Y,x=X,intercept = ii,gamma=1)
#identical(ebc1,ebc)
#which.min(ebc)

beta_ebic <- LBETA[,which.min(ebc)]

RM["MNum","EBIC"] <- sum(beta_ebic!=0)

RM["MPME","EBIC"] <- mean((Y_test-X_test%*%beta_ebic)^2)

round(RM,4)





