#A script for conducting an admixture-based genome scan

library(robust)
library(qvalue)

sv <- read.table("SV.geno", header=T) #Format: chr start end type geno1 geno2 ...
K <- read.table("SV.10.Q") #Ancestry proportions produced by ADMIXTURE
for(i in 1:ncol(K)) K[,i] <- scale(K[,i])

n <- nrow(sv)
Z <- t(sapply(1:n, function(i){
	if(i==1|i==n|i%%100==0) cat("Line",i,"/",n,"\n")
	K$y <- as.numeric(sv[i,5:ncol(sv)])
	sapply(names(K)[-11], function(x){
		as.numeric(coef(lm(formula(paste("y~",x)),K))[2])
	})
}))

S <- covRob(Z,na.action=na.omit)$cov
sv$D2 <- mahalanobis(Z,colMeans(Z,na.rm=T),S)
lambda <- median(sv$D2,na.rm=T)/qchisq(0.5,1)
sv$P <- pchisq(sv$D2/lambda,1,lower.tail=F)
sv$Q <- qvalue(sv$P)$qvalues
o <- subset(sv,Q<0.1)
o <- o[order(o$Q),]
write.table(data.frame(v1=o$chr,v2=o$start,v3=o$end),file="out.txt",row.names=F,col.names=F,quote=F,sep="\t")
