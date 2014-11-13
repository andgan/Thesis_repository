#----------------------------------------------
# Filename: Read_res.R
# Study: CCNCC
# Author: Andrea Ganna
# Date: 08OCT2010
# Updated: 16SEP2011 - Final review, picture are modified
# Purpose: Read Results from All_paper and classic_pm
# Note: This pgm needs to be run for each BIO and each sampling ratio to obtain  the results
#-----------------------------------------------
# Data used: 1X_bioX.Rdata
# Data created: 
#-----------------------------------------------
# OP: R 2.12.1
#-----------------------------------------------*/


#setwd("E:/Phd_KI/CCNCC/Results")
setwd("/Users/AndreaGanna/Documents/Work/Phd_KI/CCNCC/Results")

#############
### TABLE 2.#
#############


# Empirical efficency 
F11_beta_M3 <- NULL
for (i in 1:9){
F11_beta_m <- round(exp(mean(save[[1]][,i])),2)
F11_eff <- round(((sd(save[[1]][,1])/sd(save[[1]][,i]))^2)*100,0)
temp <- paste(F11_beta_m," ","(", F11_eff,"%",")",sep="")
F11_beta_M3 <- c(F11_beta_M3,temp)}
dim(F11_beta_M3)  <- c(1,9)
colnames(F11_beta_M3) <- names(save[[1]][1:9])

#Size
F11_size_m <-round(sapply(save[[3]],mean))


### INDIVIDUAL RISK

# Mean

mF11_Pr <- cbind(aggregate(save[[4]],list(TWINNR=save[[4]]$TWINNR),mean))
mF11_Pr_wc <- cbind(aggregate(save[[5]],list(TWINNR=save[[5]]$TWINNR),mean))
mF11_ncc <- cbind(aggregate(save[[6]],list(TWINNR=save[[6]]$TWINNR),mean))
mF11_ncc_wc <- cbind(aggregate(save[[7]],list(TWINNR=save[[7]]$TWINNR),mean))
mF11_BII <- cbind(aggregate(save[[8]],list(TWINNR=save[[8]]$TWINNR),mean))
mF11_BII_wc <- cbind(aggregate(save[[9]],list(TWINNR=save[[9]]$TWINNR),mean))
mF11_ncc_s <- cbind(aggregate(save[[10]],list(TWINNR=save[[10]]$TWINNR),mean))
mF11_ncc_s_wc <- cbind(aggregate(save[[11]],list(TWINNR=save[[11]]$TWINNR),mean))

##############################
## TABLE 3. & SUPPL TABLE 2.##
##############################

# Differences and percentiles

diff_Pr1 <- round(median(abs(mF11_Pr[,3]-mF11_Pr[,4])),2)
diff_Pr1q <- round(as.numeric(quantile(abs(mF11_Pr[,3]-mF11_Pr[,4]),probs=c(0.05,0.95))),2)
Diff_Pr1 <- paste(diff_Pr1," ","(", diff_Pr1q[1],";",diff_Pr1q[2],")",sep="")

diff_Pr1_wc <- round(median(abs(mF11_Pr_wc[,3]-mF11_Pr_wc[,4])),2)
diff_Pr1q_wc <- round(as.numeric(quantile(abs(mF11_Pr_wc[,3]-mF11_Pr_wc[,4]),probs=c(0.05,0.95))),2)
Diff_Pr1_wc <- paste(diff_Pr1_wc," ","(", diff_Pr1q_wc[1],";",diff_Pr1q_wc[2],")",sep="")

diff_ncc1 <- round(median(abs(mF11_ncc[,3]-mF11_ncc[,4])),2)
diff_ncc1q <- round(as.numeric(quantile(abs(mF11_ncc[,3]-mF11_ncc[,4]),probs=c(0.05,0.95))),2)
Diff_ncc1 <- paste(diff_ncc1," ","(", diff_ncc1q[1],";",diff_ncc1q[2],")",sep="")

diff_ncc1_wc <- round(median(abs(mF11_ncc_wc[,3]-mF11_ncc_wc[,4])),2)
diff_ncc1q_wc <- round(as.numeric(quantile(abs(mF11_ncc_wc[,3]-mF11_ncc_wc[,4]),probs=c(0.05,0.95))),2)
Diff_ncc1_wc <- paste(diff_ncc1_wc," ","(", diff_ncc1q_wc[1],";",diff_ncc1q_wc[2],")",sep="")


diff_BII1 <- round(median(abs(mF11_BII[,3]-mF11_BII[,4])),2)
diff_BII1q <- round(as.numeric(quantile(abs(mF11_BII[,3]-mF11_BII[,4]),probs=c(0.05,0.95))),2)
Diff_BII1 <- paste(diff_BII1," ","(", diff_BII1q[1],";",diff_BII1q[2],")",sep="")

diff_BII1_wc <- round(median(abs(mF11_BII_wc[,3]-mF11_BII_wc[,4])),2)
diff_BII1q_wc <- round(as.numeric(quantile(abs(mF11_BII_wc[,3]-mF11_BII_wc[,4]),probs=c(0.05,0.95))),2)
Diff_BII1_wc <- paste(diff_BII1_wc," ","(", diff_BII1q_wc[1],";",diff_BII1q_wc[2],")",sep="")

diff_ncc_s1 <- round(median(abs(mF11_ncc_s[,3]-mF11_ncc_s[,4])),2)
diff_ncc_s1q <- round(as.numeric(quantile(abs(mF11_ncc_s[,3]-mF11_ncc_s[,4]),probs=c(0.05,0.95))),2)
Diff_ncc_s1 <- paste(diff_ncc_s1," ","(", diff_ncc_s1q[1],";",diff_ncc_s1q[2],")",sep="")


diff_ncc_s1_wc <- round(median(abs(mF11_ncc_s_wc[,3]-mF11_ncc_s_wc[,4])),2)
diff_ncc_s1q_wc <- round(as.numeric(quantile(abs(mF11_ncc_s_wc[,3]-mF11_ncc_s_wc[,4]),probs=c(0.05,0.95))),2)
Diff_ncc_s1_wc <- paste(diff_ncc_s1_wc," ","(", diff_ncc_s1q_wc[1],";",diff_ncc_s1q_wc[2],")",sep="")


### WITHOUT ABSOLUTE DIFERENCE


diff_Pr1 <- round(median(mF11_Pr[,3]-mF11_Pr[,4]),2)
diff_Pr1q <- round(as.numeric(quantile(mF11_Pr[,3]-mF11_Pr[,4],probs=c(0.05,0.95))),2)
Diff_Pr1 <- paste(diff_Pr1," ","(", diff_Pr1q[1],";",diff_Pr1q[2],")",sep="")

diff_Pr1_wc <- round(median(mF11_Pr_wc[,3]-mF11_Pr_wc[,4]),2)
diff_Pr1q_wc <- round(as.numeric(quantile(mF11_Pr_wc[,3]-mF11_Pr_wc[,4],probs=c(0.05,0.95))),2)
Diff_Pr1_wc <- paste(diff_Pr1_wc," ","(", diff_Pr1q_wc[1],";",diff_Pr1q_wc[2],")",sep="")

diff_ncc1 <- round(median(mF11_ncc[,3]-mF11_ncc[,4]),2)
diff_ncc1q <- round(as.numeric(quantile(mF11_ncc[,3]-mF11_ncc[,4],probs=c(0.05,0.95))),2)
Diff_ncc1 <- paste(diff_ncc1," ","(", diff_ncc1q[1],";",diff_ncc1q[2],")",sep="")

diff_ncc1_wc <- round(median(mF11_ncc_wc[,3]-mF11_ncc_wc[,4]),2)
diff_ncc1q_wc <- round(as.numeric(quantile(mF11_ncc_wc[,3]-mF11_ncc_wc[,4],probs=c(0.05,0.95))),2)
Diff_ncc1_wc <- paste(diff_ncc1_wc," ","(", diff_ncc1q_wc[1],";",diff_ncc1q_wc[2],")",sep="")


diff_BII1 <- round(median(mF11_BII[,3]-mF11_BII[,4]),2)
diff_BII1q <- round(as.numeric(quantile(mF11_BII[,3]-mF11_BII[,4],probs=c(0.05,0.95))),2)
Diff_BII1 <- paste(diff_BII1," ","(", diff_BII1q[1],";",diff_BII1q[2],")",sep="")

diff_BII1_wc <- round(median(mF11_BII_wc[,3]-mF11_BII_wc[,4]),2)
diff_BII1q_wc <- round(as.numeric(quantile(mF11_BII_wc[,3]-mF11_BII_wc[,4],probs=c(0.05,0.95))),2)
Diff_BII1_wc <- paste(diff_BII1_wc," ","(", diff_BII1q_wc[1],";",diff_BII1q_wc[2],")",sep="")

diff_ncc_s1 <- round(median(mF11_ncc_s[,3]-mF11_ncc_s[,4]),2)
diff_ncc_s1q <- round(as.numeric(quantile(mF11_ncc_s[,3]-mF11_ncc_s[,4],probs=c(0.05,0.95))),2)
Diff_ncc_s1 <- paste(diff_ncc_s1," ","(", diff_ncc_s1q[1],";",diff_ncc_s1q[2],")",sep="")


diff_ncc_s1_wc <- round(median(mF11_ncc_s_wc[,3]-mF11_ncc_s_wc[,4]),2)
diff_ncc_s1q_wc <- round(as.numeric(quantile(mF11_ncc_s_wc[,3]-mF11_ncc_s_wc[,4],probs=c(0.05,0.95))),2)
Diff_ncc_s1_wc <- paste(diff_ncc_s1_wc," ","(", diff_ncc_s1q_wc[1],";",diff_ncc_s1q_wc[2],")",sep="")



##############
## TABLE 4. ##
##############

### PREDICTION MEASURES

### GOF ###

### OLD METHOD ###
chis_4df <- (4/4)^(1/3)
sd_4df <- ((1/(3^2))*(2/4))^(1/2)
m_4df <- 1-((1/(3^2))*(2/4))

chis_wc <- (4.75846827090/4)^(1/3)
m_wc <- 1-((1/(3^2))*(2/4))



chis_or <- as.numeric(save[[12]][,1])
chis_or_norm <- (chis_or/4)^(1/3)
chis_M<- aggregate(chis_or,list(save[[12]][,2]),mean, na.rm=T)


chis_M_norm <- (chis_M$x/4)^(1/3)
rbind(2*pnorm(-abs((chis_M_norm-chis_wc)/sd_4df)),chis_M$Group.1)


rej <- ifelse(2*pnorm(-abs((chis_or_norm-chis_wc)/sd_4df))<0.05,1,0)
a <- aggregate(rej,list(save[[12]][,2]),sum, na.rm=T)
b <- aggregate(save[[12]][,1][!is.na(save[[12]][,1])],list(save[[12]][,2][!is.na(save[[12]][,1])]),length)
a$x/b$x

### NEW EASY METHOD ###
num <- as.numeric(save[[12]][,1])
names <- save[[12]][,2]
library(VIM)
aggregate(num,list(names),countNA)
res <- aggregate(num,list(names),mean, na.rm=T)

1-pchisq(res$x,df=4)



### NRI ###

# Average
mF11_NRI <- as.matrix(aggregate(save[[14]][,1],list(name=save[[14]][,3]),mean))

# Empirical
s2F11_NRI <- as.matrix(aggregate(save[[14]][,1],list(name=save[[14]][,3]),sd))

# Empirical efficency 
F11_NRI_M3 <- NULL
for (i in 1:5){
F11_beta_m <- round(as.numeric(mF11_NRI[i,2]),1)
F11_eff <- round((((as.numeric(s2F11_NRI[1,2])/as.numeric(s2F11_NRI[i,2]))^2))*100,0)
temp <- paste(F11_beta_m," ","(", F11_eff,"%",")",sep="")
F11_NRI_M3 <- c(F11_NRI_M3,temp)}
dim(F11_NRI_M3)  <- c(1,5)
colnames(F11_NRI_M3) <- mF11_NRI[,1]


### CIND ###

#Average
mF11_CIND <- aggregate(as.numeric(save[[16]][,1]),list(name=save[[16]][,2]),mean,na.rm=T)

# Empirical var
s2F11_CIND <- as.matrix(aggregate(save[[16]][,1],list(name=save[[16]][,2]),sd))

#Empirical efficency
F11_CIND_M <- NULL
for (i in 1:10){
if (mF11_CIND$name[i] %in% c("CIND_wc_wb","CIND_Pr_wb","CIND_ncc_wb","CIND_BII_wb","CIND_ncc_s_wb"))
{F11_beta_m <- round(as.numeric(mF11_CIND[i,2]),3)
F11_eff <- round((((as.numeric(s2F11_CIND[1,2])/as.numeric(s2F11_CIND[i,2]))^2))*100,0)
temp <- paste(F11_beta_m," ","(", F11_eff,"%",")",sep="")
F11_CIND_M <- c(F11_CIND_M,temp)}}
dim(F11_CIND_M)  <- c(1,5)
colnames(F11_CIND_M) <- c("CIND_wc_wb","CIND_Pr_wb","CIND_ncc_wb","CIND_BII_wb","CIND_ncc_s_wb")


##############
## TABLE 5. ##
##############


#Average
mF11_CIND <- aggregate(as.numeric(save[[16]][,1]),list(name=save[[16]][,2]),mean,na.rm=T)

# Average of differences
diff <- NULL
for (i in seq(from=1, to=nrow(save[[16]]),by=2)){
	temp <- save[[16]][i,1] - save[[16]][i+1,1]
	diff <- c(diff,temp)}

## CHANGE to 200 for old version
names <- rep(c("CIND_wc","CIND_Pr","CIND_ncc","CIND_BII","CIND_ncc_s"),2000)
dF11_CIND <- aggregate(diff,list(names),mean,na.rm=T)

#Variability of differences
 sddF11_CIND <- aggregate(diff,list(names),sd,na.rm=T)

#Test if significant
pF11_CIND <- 2*pnorm(-abs(dF11_CIND[,2]/sddF11_CIND[,2]))
names(pF11_CIND) <- c("CIND_BII","CIND_ncc","CIND_ncc_s","CIND_Pr","CIND_wc")




### GRAPHS


mF11_Pr$c<-ifelse(mF11_Pr[,5]<0.1,21,4)
mF11_ncc$c<-ifelse(mF11_ncc[,5]<0.1,21,4)
mF11_BII$c<-ifelse(mF11_BII[,5]<0.1,21,4)
mF11_ncc_s$c<-ifelse(mF11_ncc_s[,5]<0.1,21,4)

mF11_Pr$co<-ifelse(mF11_Pr[,5]<0.1,"darkgray","black")
mF11_ncc$co<-ifelse(mF11_ncc[,5]<0.1,"darkgrey","black")
mF11_BII$co<-ifelse(mF11_BII[,5]<0.1,"darkgrey","black")
mF11_ncc_s$co<-ifelse(mF11_ncc_s[,5]<0.1,"darkgrey","black")

jpeg("Figure2_A.jpg", width=1200, height=1200, pointsize = 25, quality=100)

par(oma=c(0.5,0.5,0.5,0.5))
plot (mF11_Pr[,3], mF11_Pr[,4], pch=mF11_Pr$c,cex=1.2, col=mF11_Pr$co, log="xy", xlim=c(0.05,70), ylim=c(0.05,70), xlab="Realizations - Log Individual Risk", ylab="Case-cohort - Log Individual Risk", xaxt="n" , yaxt="n", axes=F)
mtext("A)",NORTH<-3, outer=TRUE, adj=0, padj=5)
axis(1,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)
axis(2,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)

abline(0,1)

dev.off()

jpeg("Figure2_B.jpg", width=1200, height=1200, pointsize = 25, quality=100)

par(oma=c(0.5,0.5,0.5,0.5))
plot (mF11_ncc[,3], mF11_ncc[,4], pch=mF11_ncc$c,cex=1.2, col=mF11_ncc$co, log="xy", xlim=c(0.05,70), ylim=c(0.05,70), xlab="Realizations - Log Individual Risk", ylab="Case-control - Log Individual Risk", xaxt="n" , yaxt="n", axes=F)
mtext("B)",NORTH<-3, outer=TRUE, adj=0, padj=5)
axis(1,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)
axis(2,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)

abline(0,1)

dev.off()

jpeg("Figure2_C.jpg", width=1200, height=1200, pointsize = 25, quality=100)

par(oma=c(0.5,0.5,0.5,0.5))
plot (mF11_BII[,3], mF11_BII[,4], pch=mF11_BII$c,cex=1.2, col=mF11_BII$co, log="xy", xlim=c(0.05,70), ylim=c(0.05,70), xlab="Realizations - Log Individual Risk", ylab="Case-cohort - Log Individual Risk", xaxt="n" , yaxt="n", axes=F)
mtext("C)",NORTH<-3, outer=TRUE, adj=0, padj=5)
axis(1,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)
axis(2,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)

abline(0,1)

dev.off()


jpeg("Figure2_D.jpg", width=1200, height=1200, pointsize = 25, quality=100)

par(oma=c(0.5,0.5,0.5,0.5))
plot (mF11_ncc_s[,3], mF11_ncc_s[,4], pch=mF11_ncc_s$c,cex=1.2, col=mF11_ncc_s$co, log="xy", xlim=c(0.05,70), ylim=c(0.05,70), xlab="Realizations - Log Individual Risk", ylab="Case-control - Log Individual Risk", xaxt="n" , yaxt="n", axes=F)
mtext("D)",NORTH<-3, outer=TRUE, adj=0, padj=5)
axis(1,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)
axis(2,at=c(0.001,1,10,20,40,70), las=1, label=c(0.001,1,10,20,40,70), col='black', lwd=1)

abline(0,1)

dev.off()




