#A script for testing the association between SV genotypes and gene expression
#Tuomas Hämälä, April 2021

library(qvalue)

#Common variants

df <- read.table("common.exon.exp.txt", header=T) #File produced by SV2exp.c
pc <- read.table("PC.txt", header=T) #10 PCs for controlling population structure. Format: id PC1 PC2 ...
df <- merge(df, pc, by="id")

g <- unique(df$n)
o <- data.frame(matrix(nrow=length(g),ncol=8))
colnames(o) <- c("gene","type","dist","f","b","r2","p","q")
j <- 1
for(i in 1:length(g)){
	if(i==1|i==length(g)|i%%100==0) cat("Gene",i,"/",length(g),"\n")
	temp <- subset(df,n==g[i])
	tryCatch({
		o$gene[j] <- as.character(temp$gene[1])
		o$type[j] <- as.character(temp$type[1])
		o$dist[j] <- temp$dist[1]
		m0 <- lm(scale(rna)~scale(dna)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=temp)
		m1 <- lm(scale(rna)~scale(dna)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+geno,data=temp)
		o$b[j] <- coef(m1)[13]
		o$r2[j] <- 1-deviance(m1)/deviance(m0)
		o$p[j] <- unlist(anova(m0,m1)[6])[2]
		o$f[j] <- sum(temp$geno)/(nrow(temp)*2)
		if(o$f[j] > 0.5) o$f[j] <- 1 - o$f[j]
		j <- j + 1
	},error=function(e){cat("Error in line",i,":",conditionMessage(e),"\n")})
}
o <- na.omit(o)
o$q <- qvalue(o$p)$qvalues
o <- o[order(o$q),]
nrow(subset(o,q<0.1))
write.table(o,file="common.out.txt",row.names=F,col.names=F,quote=F,sep="\t")

#Rare variants

df <- read.table("rare.exon.exp.txt", header=T) #File produced by SV2exp.c
pc <- read.table("PC.txt", header=T) #10 PCs for controlling population structure. Format: id PC1 PC2 ...
df <- merge(df, pc, by="id")

g <- unique(df$n)
o <- data.frame(matrix(nrow=length(g),ncol=6))
colnames(o) <- c("gene","type","dist","f","Z","p")
j <- 1
for(i in 1:length(g)){
	if(i==1|i==length(g)|i%%100==0) cat("Gene",i,"/",length(g),"\n")
	temp <- subset(df,n==g[i])
	temp$rna <- resid(lm(rna~dna+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=temp))
	median <- median(subset(temp,geno==0)$rna)
	mad <- mad(subset(temp,geno==0)$rna)
	temp$rna <- (temp$rna-median)/mad
	tryCatch({
		o$gene[j] <- as.character(temp$gene[1])
		o$type[j] <- as.character(temp$type[1])
		o$dist[j] <- temp$dist[1]
		o$Z[j] <- median(subset(temp,geno==1|geno==2)$rna)
		o$p[j] <- 2*pnorm(abs(o$Z[j]),lower.tail=F)
		o$f[j] <- sum(temp$geno)/(nrow(temp)*2)
		if(o$f[j] > 0.5) o$f[j] <- 1 - o$f[j]
		j <- j + 1
	},error=function(e){cat("Error in line",i,":",conditionMessage(e),"\n")})
}
o <- na.omit(o)
o$q <- qvalue(o$p)$qvalues
o <- o[order(o$q),]
nrow(subset(o,q<0.1))
write.table(o,file="rare.out.txt",row.names=F,col.names=F,quote=F,sep="\t")
