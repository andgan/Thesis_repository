#----------------------------------------------
# Filename: Bottst_WC.R
# Study: CCNCC
# Author: Andrea Ganna
# Date: 08OCT2010
# Updated: 04MAY2011 - one part was missed and now 2000 simulations
#          11AUG2011 - Now introduced the confidence intervals
# Purpose: Bottstrapped HR and SE for the realizations
# Note: Table 1 in the paper
#-----------------------------------------------
# Data used: wc.csv
# Data created: 
#-----------------------------------------------
# OP: R 2.12.1, survival 2.35-8, survey 3.22-4, gdata 2.8.0, ipred 0.8-8, aod 1.2
#-----------------------------------------------*/



nrep <- 2000


### SETTING WORKING DIRECTORY ###

#setwd("E:Phd_KI/CCNCC")
setwd("/Users/AndreaGanna/Documents/Work/Phd_KI/CCNCC")



### PACKAGES ###
#Charge packages
library(survival)




### READ FILE ###
wc_all <- read.csv("Data/wc.csv")

### CREATING SIMULATED BIOMARKER ###


# Highly associated and not correlated biomarker
set.seed(123)
temp <- ifelse(wc_all$CVD==0,rnorm(nrow(wc_all[wc_all$CVD==0,]),mean=3,sd=6),rnorm(nrow(wc_all[wc_all$CVD==1,]),mean=6,sd=6))

# Standardize
wc_all$bio1 <- temp 


# Less associated and correlated biomarker

# Generate marker
set.seed(123)
temp <- ifelse(wc_all$CVD==0,rnorm(nrow(wc_all[wc_all$CVD==0,]),mean=3,sd=6),rnorm(nrow(wc_all[wc_all$CVD==1,]),mean=4.0,sd=6))

#Set correlation with age (0.51 with age)
temp <- 0.3*wc_all$age+temp*(1-0.3)^0.5

# Set correlation with other variables (0.45 with sbp, 0.45 with antyhyp)
X <- c(wc_all$sbp,wc_all$antihyp,temp)
dim(X) <- c(nrow(wc_all),3)
rho <- 0.1
M <- array(rho,dim=c(3,3))
diag(M) <- 1
cF <- chol(M)
Y <- X %*% cF
temp2 <- Y[,3]

#Standardize
wc_all$bio2 <- temp2


################################################
### BASELINE CHARACTERISTICS AND ASSOCIATION ###
################################################

# Risk factors distribution
riskfmen<-aggregate(data.frame(wc_all$age,wc_all$TC,wc_all$sbp,wc_all$bio1,wc_all$bio2),by=list(wc_all$SEX),mean)
riskfsd<-aggregate(data.frame(wc_all$age,wc_all$TC,wc_all$sbp,wc_all$bio1,wc_all$bio2),by=list(wc_all$SEX),sd)
riskcon <- xtabs(cbind(smoke,antihyp,diabetess)~SEX, data=wc_all)


### STANDARDIZATIONS ###
wc_all$bio2 <- (wc_all$bio2-mean(wc_all$bio2))/sd(wc_all$bio2)
wc_all$bio1 <- (wc_all$bio1-mean(wc_all$bio1))/sd(wc_all$bio1)

#Association
ass <-coxph(Surv(start,stop,CVD)~age + TC  + SEX + HDL + smoke + antihyp +sbp + diabetess, data=wc_all)
assbio1 <-coxph(Surv(start,stop,CVD)~age + TC  + SEX + smoke + antihyp +sbp + diabetess + bio1, data=wc_all)
assbio2 <-coxph(Surv(start,stop,CVD)~age + TC  + SEX + smoke + antihyp +sbp + diabetess + bio2, data=wc_all)



ccfun = function (n,m,t) {

#####################################################
### WHOLE COHORT IS BOOSTRAPPED WITH REPLACEMENT ####
#####################################################

#Set seed to have reproducible results
set.seed(123+k)

wc <- wc_all[sample(nrow(wc_all),replace=T),]

#We rename TWINNR so no problem with dupicated
wc$TWINNR_or <- wc$TWINNR
wc$TWINNR <- 1:nrow(wc)


#We break ties at random
dupl <-  duplicated(wc$stop)
random<-rnorm(nrow(wc),0,1)/50
wc$stop <- ifelse(dupl,wc$stop+random,wc$stop)

########################
### 1. WHOLE COHORT ####
########################

#Set 3 year time
t <- 1095


#Set data
wc_s<-wc

### a. ASSOCIATION ###

#With Biomarker
WC_bio1 <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio1 + diabetess, data=wc_s)
seWC_bio1<- sqrt(diag(WC_bio1$var))
names(seWC_bio1) <- names(WC_bio1$coef)
coWC_bio1<- WC_bio1$coef


WC_bio2 <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio2 + diabetess, data=wc_s)
seWC_bio2<- sqrt(diag(WC_bio2$var))
names(seWC_bio2) <- names(WC_bio2$coef)
coWC_bio2<- WC_bio2$coef


WC_hdl <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + HDL + diabetess, data=wc_s)
seWC_hdl<- sqrt(diag(WC_hdl$var))
names(seWC_hdl) <- names(WC_hdl$coef)
coWC_hdl<- WC_hdl$coef

return(list(seWC_bio1,coWC_bio1,seWC_bio2,coWC_bio2,seWC_hdl,coWC_hdl))
}


F11_seWC_bio1 <- NULL
F11_coWC_bio1 <- NULL
F11_seWC_bio2 <- NULL
F11_coWC_bio2 <- NULL
F11_seWC_hdl <- NULL
F11_coWC_hdl <- NULL


for (k in 1:nrep) { F11<- ccfun(229,2,1095)
	                F11_seWC_bio1 <-rbind(F11_seWC_bio1, F11[[1]])
	                F11_coWC_bio1 <-rbind(F11_coWC_bio1, F11[[2]])
	                F11_seWC_bio2 <-rbind(F11_seWC_bio2, F11[[3]])
	                F11_coWC_bio2 <-rbind(F11_coWC_bio2, F11[[4]])
	                F11_seWC_hdl <-rbind(F11_seWC_hdl, F11[[5]])
	                F11_coWC_hdl <-rbind(F11_coWC_hdl, F11[[6]])}

# BIO1	                
round(apply(F11_seWC_bio1,2,mean),2)	
round(exp(apply(F11_coWC_bio1,2,mean)),2)	


up <- round(exp(apply(F11_coWC_bio1,2,mean)+1.959964*apply(F11_seWC_bio1,2,mean)),2)
down <- round(exp(apply(F11_coWC_bio1,2,mean)-1.959964*apply(F11_seWC_bio1,2,mean)),2)
                
# P- value for smoke that looks dfferent
2*pnorm(-abs(apply(F11_coWC_bio1,2,mean)[4]/apply(F11_seWC_bio1,2,mean)[4]))

#BIO2
round(apply(F11_seWC_bio2,2,mean),2)	
round(exp(apply(F11_coWC_bio2,2,mean)),2)

up <- round(exp(apply(F11_coWC_bio2,2,mean)+1.959964*apply(F11_seWC_bio2,2,mean)),2)
down <- round(exp(apply(F11_coWC_bio2,2,mean)-1.959964*apply(F11_seWC_bio2,2,mean)),2)
	                

#HDL
round(apply(F11_seWC_hdl,2,mean),2)	
round(exp(apply(F11_coWC_hdl,2,mean)),2)

up <- round(exp(apply(F11_coWC_hdl,2,mean)+1.959964*apply(F11_seWC_hdl,2,mean)),2)
down <- round(exp(apply(F11_coWC_hdl,2,mean)-1.959964*apply(F11_seWC_hdl,2,mean)),2)
	    
	                