#
# Filename:      GOF.R
# Content:       Find the average value of GOF in the REAL POPULATION.
#                DESIGNS:
#                 1)Whole cohort 2)Unstratified case cohort 3)Unmathced nested case control
#                 4)Stratified cohort - 4 strata 5)Individual matching nested case control 
#                FOR EACH DESIGN WE CALCULATED:
#                 a) Association b)Individual risk c)GOF 
#
# Author:        Andrea Ganna
# Created:       20 jan 2011
# Edited:        20 jan 2011 
#


nrep <- 5


### SETTING WORKING DIRECTORY ###

setwd("/Users/AndreaGanna/Documents/Work/Phd_KI/CCNCC")
#setwd("E:/Phd_KI/CCNCC")


### PACKAGES ###
#Charge packages
library(foreign)
library(survival)
library(survey)
library(gdata)
library(ipred)
library(aod)


### READ FILE ###
wc_all <- read.csv("Data/wc.csv")


### LOAD FUNCTIONS ###


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



### c.PREDICTION MEASURES ###


#GOF
wc_all$cutRt <- as.factor(cut(log(Rt_go_wb), quantile(log(Rt_go_wb),(0:7)/7), labels=1:7, include.lowest=T))
GOLD2 <- coxph(Surv(start,stop,CVD)~age + TC  + SEX  +  smoke + antihyp +sbp + bio + diabetess + cutRt, data=wc_all)

GOF<- try(wald.test(vcov(GOLD2),coef(GOLD2), Terms=9:14),silent=T)
GOF_go_wb <- try(cbind(as.vector(GOF$result$chi2)[1],"GOF_go_wb"),silent=T)
colnames(GOF_go_wb) <- c("GOF","name")

GOF_go_wb_ll <- cbind(summary(GOLD2)$logtest[1]-summary(GOLD)$logtest[1],"GOF_go_wb_ll")


 
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



### b. INDIVIDUAL RISK ###

#XB
cf <- WC$coefficient
expXB_wc <- exp(cf[1]*wc_s$age + cf[2]*wc_s$TC + cf[3]*wc_s$SEX  + cf[4]*wc_s$smoke + cf[5]*wc_s$antihyp + cf[6]*wc_s$sbp + cf[7]*wc_s$bio + cf[8]*+wc_s$diabetess)

#Baseline Hazard
H0_wc <- basehaz(WC, center=FALSE)


##Find the t-years H0	
St_wc <- exp(-H0_wc$hazard[which(abs(H0_wc$time-t)==min(abs(H0_wc$time-t)))])
	
#Individual risk
Rt_wc <- (1-St_wc^expXB_wc)*100

Rt_wc <- data.frame(cbind(Rt_wc,wc_s$TWINNR_or))
colnames(Rt_wc) <- c("Risk at t-years","TWINNR")

#For prediction measure
Rt_wc_wb <- (1-St_wc^expXB_wc)


####


#GOF
wc_s$cutRt <- as.factor(cut(log(Rt_wc_wb), quantile(log(Rt_wc_wb),(0:7)/7), labels=1:7, include.lowest=T))
tag <- 0
withCallingHandlers(
WC2 <- coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio + diabetess + cutRt, data=wc_s),
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(u <- vcov(WC2),coef(WC2), Terms=9:14)
GOF_wc_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_wc_wb")
colnames(GOF_wc_wb) <- c("GOF","name")} else {GOF_wc_wb <- c(NA,"GOF_wc_wb")}

GOF_wc_wb_ll <- cbind(summary(WC2)$logtest[1]-summary(WC)$logtest[1],"GOF_wc_wb_ll")



### NOW START WITH SAMPLING DESIGN ###


####################################
### 2. UNSTRATIFIED CASE COHORT ####
#################################### 

#Set data
wc_s<-wc

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

#Prentice
PR <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess + cluster(as.factor(keys)), data=cc_s)

sePR <- sqrt(diag(PR$var)[7])
coPR <- PR$coef[7]


### b. INDIVIDUAL RISK ###

#XB
cf <- PR$coefficient
cc_s$expXB_Pr <- exp(cf[1]*cc_s$age + cf[2]*cc_s$TC + cf[3]*cc_s$SEX + cf[4]*cc_s$smoke + cf[5]*cc_s$antihyp + cf[6]*cc_s$sbp + cf[7]*cc_s$bio + cf[8]*cc_s$diabetess )

#Creating failures time 
failure <- cc_s[which(cc_s$CVD_cc==1),]
failure <- failure[order(failure$stop),]
  
#Define risk sets
sexpXB_Pr <- NULL
for (a in 1:dim(failure)[1]) {
	B <- cc_s[which(cc_s$start_cc < failure$stop[a] & failure$stop[a] <= cc_s$stop_cc),] 
	sexpXB_Pr <- c( sexpXB_Pr,sum(B$expXB_Pr))}
   
#Baseline Hazard
H0_Pr<- data.frame(cbind(failure$stop,(cumsum(prob*(1/sexpXB_Pr)))))
colnames(H0_Pr) <- c("time","hazard")

#Find the t-years	
St_Pr <- exp(-H0_Pr$hazard[which(abs(H0_Pr$time-t)==min(abs(H0_Pr$time-t)))])
	
#Individual Risk
Rt_Pr <- (1-St_Pr^cc_s$expXB_Pr)*100

Rt_Pr <- data.frame(cbind(Rt_Pr,cc_s$TWINNR,cc_s$CVD_cc))
colnames(Rt_Pr) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_Pr_wb <- (1-St_Pr^cc_s$expXB_Pr)


#GOF
cc_s$cutRt <- as.factor(cut(log(Rt_Pr_wb), quantile(log(Rt_Pr_wb),(0:7)/7), labels=1:7, include.lowest=T))
tag <- 0
withCallingHandlers(
PR2 <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess +   cutRt + cluster(as.factor(keys)), data=cc_s),
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(vcov(PR2),coef(PR2), Terms=9:14)
GOF_Pr_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_Pr_wb")
colnames(GOF_Pr_wb) <- c("GOF","name")} else {GOF_Pr_wb <- c(NA,"GOF_Pr_wb")}




#########################################
### 3. UNMATCHED NESTED CASE CONTROL ####
######################################### 

#Set data
wc_s<-wc
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

### b.INDIVIDUAL RISK ###

#XB  
cf <- NCC$coefficient
expXB_ncc <- exp(cf[1]*ncc$age + cf[2]*ncc$TC + cf[3]*ncc$SEX  + cf[4]*ncc$smoke + cf[5]*ncc$antihyp + cf[6]*ncc$sbp + cf[7]*ncc$bio +cf[8]*ncc$diabetess)

#Baseline Hazard
H0_ncc <- basehaz(NCC, center=FALSE)

#Find the t-years 
St_ncc <- exp(-H0_ncc$hazard[which(abs(H0_ncc$time-t)==min(abs(H0_ncc$time-t)))])
	
#Individual risk	
Rt_ncc <- (1-St_ncc^expXB_ncc)*100
Rt_ncc <- data.frame(cbind(Rt_ncc,ncc$TWINNR,ncc$cc))
colnames(Rt_ncc) <- c("Risk at t-years","TWINNR")

# For prediction measures
Rt_ncc_wb <- (1-St_ncc^expXB_ncc)

#GOF
ncc$cutRt <- as.factor(cut(log(Rt_ncc_wb), quantile(log(Rt_ncc_wb),(0:7)/7), labels=1:7, include.lowest=T))
tag <- 0
withCallingHandlers(
NCC2 <- coxph(Surv(fakentry,rstime,cc)~age+TC+ SEX + smoke + antihyp + sbp + bio + diabetess + cutRt + offset(logw), data=ncc), 
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(vcov(NCC2),coef(NCC2), Terms=9:14)
GOF_ncc_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_ncc_wb")
colnames(GOF_ncc_wb) <- c("GOF","name")} else {GOF_ncc_wb <- c(NA,"GOF_ncc_wb")}

GOF_ncc_wb_ll <- cbind(summary(NCC2)$logtest[1]-summary(NCC)$logtest[1],"GOF_ncc_wb_ll")




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


### b.INDIVIDUAL RISK ###

#XB
cf <- BII_s4$coefficient
expXB_BII <- exp(cf[1]*cc$age + cf[2]*cc$TC + cf[3]*cc$SEX  + cf[4]*cc$smoke + cf[5]*cc$antihyp + cf[6]*cc$sbp + cf[7]*cc$bio + cf[8]*cc$diabetess)

#Baseline Hazard
H0_BII <- basehaz(BII_s4, center=FALSE)

#Find the t-years 	
St_BII <- exp(-H0_BII$hazard[which(abs(H0_BII$time-t)==min(abs(H0_BII$time-t)))])
	
#Individual risk	
Rt_BII <- (1-St_BII^expXB_BII)*100
Rt_BII <- data.frame(cbind(Rt_BII,cc$TWINNR,cc$CVD))
colnames(Rt_BII) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_BII_wb <- (1-St_BII^expXB_BII)



#GOF
dstrat2$phase1$sample$variables$cutY <- as.factor(cut(log(Rt_BII_wb), quantile(log(Rt_BII_wb),(0:7)/7), labels=1:7, include.lowest=T))


tag <- 0
withCallingHandlers(
BII_s42 <-svycoxph(Surv(start,stop,CVD)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess + cutY,design=dstrat2), 
warning=function(w) {tag <<-1}, error=function(e) {tag<<-1})

if(tag==0) {GOF <- wald.test(vcov(BII_s42),coef(BII_s42), Terms=9:14)
GOF_BII_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_BII_wb")
colnames(GOF_BII_wb) <- c("GOF","name")} else {GOF_BII_wb<- c(NA,"GOF_BII_wb")}
 





########################################
### 5. MATCHED NESTED CASE CONTROL  ####
########################################
 
#Set data
wc_s<-wc

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



### b.INDIVIDUAL RISK ###

#XB
cf <- NCC_s$coefficient
expXB_ncc_s <- exp(cf[1]*ncc_s$TC  + cf[2]*ncc_s$smoke + cf[3]*ncc_s$antihyp + cf[4]*ncc_s$sbp + cf[5]*ncc_s$bio + cf[6]*ncc_s$diabetess)
sexpXB_ncc_s <- aggregate(expXB_ncc_s*(ncc_s$w),list(time=ncc_s$rstime),sum)


#Baseline Hazard
H0_ncc_s<- data.frame(cbind(sexpXB_ncc_s$time,(cumsum(1/sexpXB_ncc_s$x))))
colnames(H0_ncc_s) <- c("time","hazard")

#Find the t-years 
St_ncc_s <- exp(-H0_ncc_s$hazard[which(abs(H0_ncc_s$time-t)==min(abs(H0_ncc_s$time-t)))])


#Individual risk
Rt_ncc_s <- data.frame(cbind((1-(St_ncc_s^expXB_ncc_s))*100,ncc_s$TWINNR,ncc_s$cc))
colnames(Rt_ncc_s) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_ncc_s_wb <- Rt_ncc_s[,1]/100

#GOF
ncc_s$cutRt <- as.factor(cut(log(Rt_ncc_s_wb), quantile(log(Rt_ncc_s_wb),(0:7)/7), labels=1:7, include.lowest=T))
tag <- 0
withCallingHandlers(
NCC_s2 <- coxph(Surv(fakentry,rstime,cc)~TC + smoke + antihyp + sbp + bio + diabetess + cutRt + offset(logw), data=ncc_s), 
warning=function(w) {tag <<-1})

if(tag==0) {GOF <- wald.test(vcov(NCC_s2),coef(NCC_s2), Terms=7:12)
GOF_ncc_s_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_ncc_s_wb")
colnames(GOF_ncc_s_wb) <- c("GOF","name")} else {GOF_ncc_s_wb <- c(NA,"GOF_ncc_s_wb")}
      
GOF_ncc_s_wb_ll <- cbind(summary(NCC_s2)$logtest[1]-summary(NCC_s)$logtest[1],"GOF_ncc_s_wb_ll")
 


     
#################
#MERGING RESULTS#
#################

### a.ASSOCIATION ###




### c.PREDICTION MEASURES ###

GOFm <- rbind (GOF_wc_wb, GOF_wc_wb_ll, GOF_Pr_wb,GOF_ncc_wb,GOF_ncc_wb_ll,GOF_BII_wb,GOF_ncc_s_wb,GOF_ncc_s_wb_ll)


return(list(GOFm))


}


F11 <-NULL
F11_GOF <- NULL
F11_GOF2 <- NULL
F11_pred <- NULL
F11_pred2 <- NULL
F11_obs <- NULL
F11_obs2 <- NULL
F11_HL <- NULL
F11_HL2 <- NULL

F13 <-NULL
F13_GOF <- NULL
F13_GOF2 <- NULL
F13_pred <- NULL
F13_pred2 <- NULL
F13_obs <- NULL
F13_obs2 <- NULL
F13_HL <- NULL
F13_HL2 <- NULL 

for (k in 1:nrep) {   F11<- ccfun(229,2,1095)
	                  F11_GOF <- rbind(F11_GOF, F11[[1]])
                      print(k)}



save <- list(F11_GOF)
save(save,file="GOF_bio1_1_7.Rdata")


for (k in 1:nrep) {   F13<- ccfun(632,4,1095)
	                  F13_GOF <- rbind(F13_GOF, F13[[1]])
	                                      print(k)}



save <- list(F13_GOF)
save(save,file="GOF_bio1_3_7.Rdata")



	