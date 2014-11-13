#----------------------------------------------
# Filename: All_paper_reviewer.R
# Study: CCNCC
# Author: Andrea Ganna
# Date: 20APR2011
# Updated: 20APR2011
# Purpose: All main analyses for CCNCC - reviewer comments. No prediction measures are included. This program is to answer to the reviewer comment that wants to see what happens when the number of unique subjects is set equal.
# Note: The program is now set on BIO1, needs also to be run for BIO2 and HDL to have the complete analyses
#-----------------------------------------------
# Data used: wc.csv
# Data created: 11_bioX.Rdata
#               13_bioX.Rdata
#-----------------------------------------------
# OP: R 2.12.1, survival 2.35-8, survey 3.22-4, gdata 2.8.0, ipred 0.8-8, aod 1.2
#-----------------------------------------------*/


#STRUCTURE
#
#DESIGNS:
#                 1)Whole cohort 2)Unstratified CC 3)Unmathced NCC
#                 4)Stratified CC 5)Individual matching NCC
#FOR EACH DESIGN WE CALCULATED:
#                 a) Association b)Individual risk c)Prediction measures 



nrep <- 1000


### SETTING WORKING DIRECTORY ###

#setwd("E:/Phd_KI/CCNCC")
setwd("/Users/AndreaGanna/Documents/Work/Phd_KI/CCNCC")



### PACKAGES ###
#Charge packages
library(survival)
library(survey)
library(gdata)
library(ipred)
library(aod)


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



######################################################
######################################################
#####################  ANALYSIS ######################
######################################################
######################################################


wc_all$bio <- wc_all$bio1

### GOLD STANDARD ###

#Set 3 year time
t <- 1095

### a. ASSOCIATION ###

#With Biomarker
GOLD <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio + diabetess, data=wc_all)



### b. INDIVIDUAL RISK ###

#XB
cf <- GOLD$coefficient
expXB_go <- exp(cf[1]*wc_all$age + cf[2]*wc_all$TC + cf[3]*wc_all$SEX  + cf[4]*wc_all$smoke + cf[5]*wc_all$antihyp + cf[6]*wc_all$sbp + cf[7]*wc_all$bio + cf[8]*+wc_all$diabetess)

#Baseline Hazard
H0_go <- basehaz(GOLD, center=FALSE)

#Find the t-years H0	
St_go <- exp(-H0_go$hazard[which(abs(H0_go$time-t)==min(abs(H0_go$time-t)))])
	
#For prediction measure
Rt_go_wb <- (1-St_go^expXB_go)



 
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
WC <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio + diabetess, data=wc_s)
seWC<- sqrt(diag(WC$var)[7])
coWC<- WC$coef[7]



### NOW START WITH SAMPLING DESIGN ###


####################################
### 2. UNSTRATIFIED CASE COHORT ####
#################################### 
#Set data (selecting variable to speed-up calculation)
wc_s<-wc[c("TWINNR","SEX","TWINNR_or","diabetess","TC","HDL","age","sbp","bio","smoke","CVD","start","stop","antihyp")]

#Average number of people for a random sampling
#random <- (n-(table(wc_all$CVD)[2]/dim(wc_all)[1])*n)+table(wc_all$CVD)[2]

### PREPARING DATA ###
       
#Definition subcohort
subcoh<-wc_s[sample(1:nrow(wc_s),n,replace=F),]

#Indicator for cohort sample
subcoh$ind <- as.factor(1)

#Subcoh variables
subcoh <- cbind(subcoh$TWINNR, subcoh$ind)
colnames(subcoh)<- c("TWINNR","ind")

#Merge subcoh with whole cohort
wc_s <- merge(subcoh,wc_s, by="TWINNR", all.y=T)

#Setting subcoh indicator
wc_s$ind <- ifelse (is.na(wc_s$ind) == T, 0,1)

#Keep subcohort and cases
cc<-subset(wc_s, wc_s$ind==1 | wc_s$CVD==1)

#Size
L_random <- nrow(cc)

#Arrange the case-cohort dataset
e <- 0.00001
prob <- n/nrow(wc_s) 
  
  start_cc <- NULL
  stop_cc <- NULL
  CVD_cc <- NULL 
  wBR <- NULL
  wPR <- NULL
  logwSP <- NULL 
  keys <- NULL

  for (i in 1:nrow(cc)) {

      # Case outside subcohort	
      if ((cc$CVD[i]==1) & (cc$ind[i]!=1)) 
 		{CVD_cc <- c(CVD_cc, 1)
  		 wBR <- c(wBR, 1)
  		 logwSP <- c(logwSP, -100)
  		 stop_cc <- c(stop_cc, cc$stop[i])
  		 start_cc <- c(start_cc, cc$stop[i]-e)
  		 keys <- c(keys, cc$TWINNR[i])}
  		 
      # Non-case in subcohort
	  else if ((cc$CVD[i]!=1) & (cc$ind[i]==1)) 
 		{CVD_cc <- c(CVD_cc, 0) 
  		 wBR <- c(wBR, 1/prob)
  		 logwSP <- c(logwSP, 0)
  		 start_cc <- c(start_cc, cc$start[i])
  		 stop_cc <- c(stop_cc, cc$stop[i])         
         keys <- c(keys, cc$TWINNR[i])}
  
      # Case in subcohort 
	  else if ((cc$CVD[i]==1) & (cc$ind[i]==1)) 
		{CVD_cc<- c(CVD_cc, 0)
         wBR <- c(wBR, 1/prob) 
         logwSP <- c(logwSP, 0)
         start_cc <- c(start_cc, cc$start[i])
         stop_cc <- c(stop_cc, cc$stop[i]-e)
         keys <- c(keys, cc$TWINNR[i])
         CVD_cc <- c(CVD_cc, 1)
         wBR <-  c(wBR, 1)
         logwSP <-  c(logwSP, -100)
         stop_cc <- c(stop_cc, cc$stop[i])
         start_cc <- c(start_cc, cc$stop[i]-e)
         keys <- c(keys, cc$TWINNR[i])}
     } 

#Merging obtained variables
temp <- data.frame(cbind(CVD_cc, wBR, logwSP, keys, start_cc, stop_cc))

#Merging with the original
cc_s<- merge(temp, cc, by.x="keys", by.y="TWINNR", all.y=T, all.x=T)

   
### a. ASSOCIATION ###

#With Biomarker

#Breslow 
BR <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX + smoke + antihyp + sbp + bio + diabetess + cluster(as.factor(keys)),weights=wBR, data=cc_s) 

seBR <- sqrt(diag(BR$var)[7])
coBR <- BR$coef[7] 



#Prentice
PR <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess + cluster(as.factor(keys)), data=cc_s)

sePR <- sqrt(diag(PR$var)[7])
coPR <- PR$coef[7]


#Self & Prentice
SP <- coxph(Surv(start,stop,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess + offset(logwSP) + cluster(as.factor(keys)), data=cc_s) 

seSP <- sqrt(diag(SP$var)[7])
coSP <- SP$coef[7]


#########################################
### 3. UNMATCHED NESTED CASE CONTROL ####
######################################### 

#Set data (selecting variable to speed-up calculation)
wc_s<-wc[c("TWINNR","SEX","TWINNR_or","diabetess","TC","HDL","age","sbp","bio","smoke","CVD","stop","antihyp")]

### PREPARING DATA ###


#Creating failures time 
failure <- wc_s[which(wc_s$CVD==1),]
failure <- failure[order(failure$stop),]
  
ncc <- NULL
  for (a in 1:nrow(failure)) {
  	B <- wc_s[failure$stop[a] < wc_s$stop,]
  	B2 <-rbind(B[sample(nrow(B),m-1,replace=F),],failure[a,])
  	B2$cc <- ifelse(failure$TWINNR[a]==B2$TWINNR,1,0)
  	B2$w <- (nrow(B)+1)/(m)
  	B2$rstime <- B2[B2$cc==1,]$stop
  	ncc <- rbind(ncc,B2)}
  	
ncc$logw <- log(ncc$w)
ncc$fakentry <- ncc$rstime-0.0001    
  
ncc <- ncc[order(ncc$rstime,ncc$cc,ncc$stop),]

### a.ASSOCIATION ###

#With biomarker 
NCC <- coxph(Surv(fakentry,rstime,cc)~age+TC+ SEX + smoke + antihyp + sbp + bio + diabetess + offset(logw), data=ncc) 

seNCC <- sqrt(diag(NCC$var)[7])
coNCC <- NCC$coef[7]


#Size
L_random_ncc <- length(unique(ncc$TWINNR))


#############################################
### 4. STRATIFIED CASE COHORT - 4 STRATA ####
#############################################


#Set data
wc_s<-wc

# Create strata
wc_s$strata <- ifelse(wc_s$age<median(wc_s$age) & wc_s$SEX==2,1,
               ifelse(wc_s$age<median(wc_s$age) & wc_s$SEX==1,2,
               ifelse(wc_s$age>=median(wc_s$age) & wc_s$SEX==2,3,4)))



### PREPARING DATA ###

temp<-table(wc_s[wc_s$CVD==1,]$strata)
exp <- round((n/sum(temp)*temp))

subcoh <- NULL
for (a in 1:length(exp)){
	temp<- wc_s[wc_s$strata==a,]
	subcoh<- rbind(subcoh,temp[sample(1:nrow(temp),exp[a],replace=F),])}

#Indicator for cohort sample
subcoh$ind <- as.factor(1)

#Subcoh variables
subcoh <- cbind(subcoh$TWINNR, subcoh$ind)
colnames(subcoh)<- c("TWINNR","ind")

#Merge subcoh with whole cohort
wc_s <- merge(subcoh,wc_s, by="TWINNR", all.y=T)

#Setting subcoh indicator
wc_s$ind <- ifelse (is.na(wc_s$ind) == T, 0,1)

#Keep subcohort and cases
cc<-subset(wc_s, wc_s$ind==1 | wc_s$CVD==1)

#Size
L_strata4 <- nrow(cc)

          
### a.ASSOCIATION ###

#With biomarker

#Borgan I

stratsizes<-table(wc_s$strata)
subc0h <- wc_s$ind
selccoh <- with(wc_s, CVD==1|subc0h==1)
ccoh.data <- wc_s[selccoh,]

BI_s4 <- cch(Surv(start,stop,CVD) ~ age+SEX+smoke+antihyp+TC+sbp+bio+diabetess, data=ccoh.data, subcoh=~ind, id=~TWINNR, stratum=~strata, cohort.size=stratsizes, method="I.Borgan")

seBI_s4 <- sqrt(diag(BI_s4$var)[7])
coBI_s4 <- BI_s4$coef[7]
 

# Borgan II

wc_s$strata   <- ifelse(wc_s$CVD==1,0,wc_s$strata)        
subset <- wc_s$ind==1 | wc_s$CVD==1
wc_s$bio <- ifelse(subset==T,wc_s$bio,NA)

#Step 1: Predict missing covariates for all subjects (not just for those with missing covariates)
dstrat<-twophase(id=list(~TWINNR,~TWINNR),strata=list(NULL,~strata),subset=subset ,data=wc_s) 
fit.step1 <- svyglm(bio ~ age + sbp + antihyp + ApoA1, design=dstrat)
wc_s$bio_p <- predict(fit.step1,type="response",newdata=wc_s,se=F)

# Step 2: Fit an augmented dataset with risk model to get the auxiliary variable: dfbeta
calmodel<-coxph(Surv(start,stop,CVD)~age+TC+ SEX + smoke + antihyp + sbp + bio_p + diabetess,  data=wc_s)
db = resid(calmodel,"dfbeta")+1 
colnames(db)<-paste("db",1:ncol(db),sep="")
wc_s2<-cbind(wc_s,db)
dstrat2<-twophase(id=list(~TWINNR,~TWINNR),strata=list(NULL,~strata),subset=subset ,data=wc_s2)

BII_s4<-svycoxph(Surv(start,stop,CVD)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess,design=dstrat2)  

seBII_s4 <- sqrt(diag(BII_s4$var)[7])
coBII_s4 <- BII_s4$coef[7]


#Two-stage calibration

dcal<-calibrate(dstrat2,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=nrow(wc_s),colSums(db)),calfun="raking",eps=0.0001)
CA_s4<-svycoxph(Surv(start,stop,CVD)~age+TC+ SEX  + smoke + antihyp + sbp + bio+diabetess, design=dcal)

seCA_s4 <- sqrt(diag(CA_s4$var)[7])
coCA_s4 <- CA_s4$coef[7]




 
########################################
### 5. MATCHED NESTED CASE CONTROL  ####
########################################
 
#Set data (selecting variable to speed-up calculation)
wc_s<-wc[c("TWINNR","SEX","TWINNR_or","diabetess","TC","HDL","age","sbp","bio","smoke","CVD","stop","antihyp")]
 
#Round AGE
wc_s$age <- round(wc_s$age,0)

### PREPARING DATA ###

#Creating failures time (all CVD events)
failure <- wc_s[which(wc_s$CVD==1),]
failure <- failure[order(failure$stop),]
  

ncc_s <- NULL
  for (a in 1:nrow(failure)) {
  	B <- wc_s[failure$stop[a] < wc_s$stop & failure$SEX[a]==wc_s$SEX & failure$age[a]==wc_s$age,]
  	if (nrow(B)<m) {B2<-rbind(B,failure[a,])}else{B2 <-rbind(B[sample(nrow(B),m-1,replace=F),],failure[a,])}
  	B2$w2 <- (nrow(B)+1)/m
  	B2$cc <- ifelse(failure$TWINNR[a]==B2$TWINNR,1,0)
  	B2$w <- (nrow(wc_s[failure$stop[a] < wc_s$stop,])+1)/(m)
  	B2$rstime <- B2[B2$cc==1,]$stop
  	ncc_s <- rbind(ncc_s,B2)}
  	
ncc_s$logw <- log(ncc_s$w)
ncc_s$fakentry <- ncc_s$rstime-0.0001  

ncc_s <- ncc_s[order(ncc_s$rstime,ncc_s$cc,ncc_s$stop),]	



### a.ASSOCIATION ###
 
NCC_s <- coxph(Surv(fakentry,rstime,cc)~TC + smoke + antihyp + sbp + bio + diabetess, data=ncc_s) 

seNCC_s <- sqrt(diag(NCC_s$var)[5])
coNCC_s <- NCC_s$coef[5]


#Size
L_strata_ncc <- length(unique(ncc_s$TWINNR))


# We broke the matching
NCC_s_br <- glm(CVD ~age + SEX + TC + smoke + antihyp + sbp + bio + diabetess, data=ncc_s, family=binomial)
coNCC_s_br_age <- NCC_s_br$coef[2] 
coNCC_s_br_sex <- NCC_s_br$coef[3] 


      
#################
#MERGING RESULTS#
#################

### a.ASSOCIATION ###

CO <- cbind(coWC, coPR,  coSP, coBR, coNCC, coBI_s4, coBII_s4, coCA_s4,coNCC_s,coNCC_s_br_age,coNCC_s_br_sex)
SE <- cbind(seWC, sePR,  seSP, seBR, seNCC, seBI_s4, seBII_s4, seCA_s4,seNCC_s)

#NUmber of non-repeated subjects
DIM <- round(cbind(L_random, L_random_ncc,  L_strata4, L_strata_ncc),0)


return(list(CO,SE, DIM))


}


F11 <-NULL
F11_beta <-NULL
F11_se <-NULL
F11_size <-NULL
F11_Pr <- NULL
F11_Pr_wc <- NULL
F11_ncc <- NULL
F11_ncc_wc <- NULL
F11_BII <- NULL
F11_BII_wc <- NULL
F11_ncc_s <- NULL
F11_ncc_s_wc <- NULL
F11_GOF <- NULL
F11_HL <- NULL
F11_NRI <- NULL
F11_IDI <- NULL
F11_CIND <- NULL
F13 <-NULL
F13_beta <-NULL
F13_se <-NULL
F13_size <-NULL
F13_Pr <- NULL
F13_Pr_wc <- NULL
F13_ncc <- NULL
F13_ncc_wc <- NULL
F13_BII <- NULL
F13_BII_wc <- NULL
F13_ncc_s <- NULL
F13_ncc_s_wc <- NULL
F13_GOF <- NULL
F13_HL <- NULL
F13_NRI <- NULL
F13_IDI <- NULL
F13_CIND <- NULL
 

for (k in 1:nrep) {   F11<- ccfun(238,2,1095)
	                F11_beta <-data.frame(rbind(F11_beta, F11[[1]]))
	                F11_se <-data.frame(rbind(F11_se, F11[[2]]))
	                F11_size <-data.frame(rbind(F11_size, F11[[3]]))}



save <- list(F11_beta,F11_se,F11_size)
save(save,file="Results/Analysis_reviewer/unique_subjects_comparison/11_bio1.Rdata")


for (k in 1:nrep) {   F13<- ccfun(684,4,1095)
	                F13_beta <-data.frame(rbind(F13_beta, F13[[1]]))
	                F13_se <-data.frame(rbind(F13_se, F13[[2]]))
	                F13_size <-data.frame(rbind(F13_size, F13[[3]]))}



save <- list(F13_beta,F13_se,F13_size)
save(save,file="Results/Analysis_reviewer/unique_subjects_comparison/13_bio1.Rdata")

