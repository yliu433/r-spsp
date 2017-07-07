


# load GBM gene expression data ———
ge2 <- read.table("https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/AgilentG4502A_07_1",
                 header = T,sep="\t",fill = T,quote = "",stringsAsFactors = F) # 17814 by 103

ge3 <- read.table("https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/AgilentG4502A_07_2",
                 header = T,sep="\t",fill = T,quote = "",stringsAsFactors = F) # 17814 by 484

identical(sort(ge2[,1]),sort(ge3[,1])) ## TRUE

intersect(colnames(ge2),colnames(ge3))  ## "TCGA.06.0216.01"

ge3 <- ge3[match(ge2[,1],ge3[,1]),]

ge21 <- ge2[,"TCGA.06.0216.01"]
ge31 <- ge3[,"TCGA.06.0216.01"]

ge22 <- ge2[,-1];rownames(ge22) <- ge2$sample;ge22 <- t(ge22)

ge32 <- ge3[,-1];rownames(ge32) <- ge3$sample;ge32 <- t(ge32)

sam2 <- (rownames(ge22)) # 
sam2 <- sam2[sam2!=("TCGA.06.0216.01")]
sam3 <- (rownames(ge32)) # 

Xall <- ge32[,colnames(ge32)!="TP53"]

## use the second data for testing
expY_test <- ge22[,"TP53"]
Xall_test <- ge22[match(sam2,rownames(ge22)),]

identical(as.character(pheno2$sampleID[match(sam3,pheno2$sampleID)]),
          as.character(rownames(Xall)))

Yall <- log(expY)
## center the reponse and standardize the predictors

n <- length(Yall)   ### 370
Yc <- Yall-mean(Yall)
Xs <- apply(Xall,2,scale) * sqrt(n/(n-1))

sum(Xs[,1]^2)/n

corxy <- apply(Xs,2,function(x){cor(x,Yc)})

Xall1000 <- Xall[,order(abs(corxy),decreasing = TRUE)[1:1000]]

Xs1000 <- Xs[,order(abs(corxy),decreasing = TRUE)[1:1000]]

Y0_test <- log(expY_test)
X0_test <- Xall_test[,order(abs(corxy),decreasing = TRUE)[1:1000]]

identical(colnames(X0_test),colnames(Xs1000))

meanx <- apply(Xall1000,2,mean)
sdx <- apply(Xall1000,2,sd)/sqrt(n/(n-1))


Y_test <- Y0_test-mean(Yall)
X_test_tmp <- sweep(X0_test,2,meanx,"-")
X_test <- sweep(X_test_tmp,2,sdx,"/")

save(Yc,Xs1000,X_test,Y_test,file="./gbm.Rdata")






