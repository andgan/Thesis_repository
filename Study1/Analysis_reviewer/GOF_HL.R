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


nrep <- 1


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
library(Hmisc)



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

#GOF
wc_all$cutRt <- as.factor(cut(Rt_go_wb, c(0,0.05,0.20,1), labels=1:3, include.lowest=T))


res_all <- as.numeric(residuals(GOLD, type=c("martingale")))
sub_all <- wc_all$CVD-res_all

pred_all <- aggregate(sub_all,list(wc_all$cutRt),sum)

obs_all <- aggregate(wc_all$CVD,list(wc_all$cutRt),sum)

N_all <- table (wc_all$cutRt)

HL_all <- sum ((obs_all$x- pred_all$x)^2/(pred_all$x*(1-pred_all$x/N_all)))




 
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
wc_s$cutRt <- as.factor(cut(log(Rt_wc_wb), quantile(log(Rt_wc_wb),(0:5)/5), labels=1:5, include.lowest=T))


res_wc <- as.numeric(residuals(WC, type=c("martingale")))
sub_wc <- wc_s$CVD-res_wc

pred_wc <- aggregate(sub_wc,list(wc_s$cutRt),sum)

obs_wc <- aggregate(wc_s$CVD,list(wc_s$cutRt),sum)

N <- table (wc_s$cutRt)

HL_wc <- sum ((obs_wc$x- pred_wc$x)^2/(pred_wc$x*(1-pred_wc$x/N)))
HL_wc_b <- sum ((obs_wc$x- pred_wc$x)^2/(pred_wc$x))


### GOF cat
wc_s$cutRt2 <- as.factor(cut(Rt_wc_wb, c(0,0.05,0.20,1), labels=1:3, include.lowest=T))

res_wc2 <- as.numeric(residuals(WC, type=c("martingale")))
sub_wc2 <- wc_s$CVD-res_wc2

pred_wc2 <- aggregate(sub_wc2,list(wc_s$cutRt2),sum)

obs_wc2 <- aggregate(wc_s$CVD,list(wc_s$cutRt2),sum)

N2 <- table (wc_s$cutRt2)

HL_wc2 <- sum ((obs_wc2$x- pred_wc2$x)^2/(pred_wc2$x*(1-pred_wc2$x/N2)))
HL_wc2_b <- sum ((obs_wc2$x- pred_wc2$x)^2/(pred_wc2$x))



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


#Risk Excluding cases in the subcohort
Rt_Pr_wb_m <- Rt_Pr_wb[(cc_s$CVD_cc==0 & cc_s$CVD==1)==F]

#Weights for prediction measures
w <- 1/(ifelse(cc$CVD==1,1,nrow(cc[wc_s$ind==1 & wc_s$CVD==0,])/nrow(wc_s[wc_s$CVD==0,])))
 


#GOF
cc_s$cutRt <- as.factor(cut(log(Rt_Pr_wb), quantile(log(Rt_Pr_wb),(0:5)/5), labels=1:5, include.lowest=T))

res_Pr <- as.numeric(residuals(PR, type=c("martingale")))
sub_Pr <- cc_s$CVD_cc-res_Pr

pred_Pr <- aggregate(sub_Pr,list(cc_s$cutRt),sum)

obs_Pr <- aggregate(cc_s$CVD_cc,list(cc_s$cutRt),sum)

N_Pr <- table(cc_s$cutRt)

HL_Pr <- sum ((obs_Pr$x- pred_Pr$x)^2/(pred_Pr$x))



cc_s$cutRt2 <- as.factor(cut(Rt_Pr_wb, c(0,0.5,0.20,1), labels=1:3, include.lowest=T))

res_Pr2 <- as.numeric(residuals(PR, type=c("martingale")))
sub_Pr2 <- cc_s$CVD_cc-res_Pr2

pred_Pr2 <- aggregate(sub_Pr,list(cc_s$cutRt2),sum)

obs_Pr2 <- aggregate(cc_s$CVD_cc,list(cc_s$cutRt2),sum)

if (sum(as.numeric(pred_Pr2$Group.1))!=6){
	pred_Pr2 <- rbind(pred_Pr2,c(3,0))
	obs_Pr2 <- rbind(obs_Pr2,c(3,0))}



N_Pr2 <- table(cc_s$cutRt2)

HL_Pr2 <- sum ((obs_Pr2$x- pred_Pr2$x)^2/(pred_Pr2$x))




#GOF
cc_s$cutRt3 <- as.factor(cut(log(Rt_Pr_wb), wtd.quantile(log(Rt_Pr_wb_m),weights=w,(0:5)/5), labels=1:5, include.lowest=T))
cc_s$cutRt3 <- ifelse(log(Rt_Pr_wb)>wtd.quantile(log(Rt_Pr_wb_m),weights=w,(0:5)/5)[6],5,
               ifelse(log(Rt_Pr_wb)<wtd.quantile(log(Rt_Pr_wb_m),weights=w,(0:5)/5)[1],1,cc_s$cutRt3))

res_Pr3 <- as.numeric(residuals(PR, type=c("martingale")))
sub_Pr3 <- cc_s$CVD_cc-res_Pr3

pred_Pr3 <- aggregate(sub_Pr3,list(cc_s$cutRt3),sum)

obs_Pr3 <- aggregate(cc_s$CVD_cc,list(cc_s$cutRt3),sum)

N_Pr3 <- table(cc_s$cutRt3)

HL_Pr3 <- sum ((obs_Pr3$x- pred_Pr3$x)^2/(pred_Pr3$x))






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


temp<- ncc[c("stop","cc")]
w=rep(1,nrow(ncc))
for(i in 1:nrow(temp)){
	if(temp$cc[i]==0) {
		time <- temp$stop[temp$cc==1 & temp$stop <= temp$stop[i]]
		prob.pvs <- rep(NA,length(time))
		for(j in 1:length(time)) {
			prob.pvs[j] <- (m-1)/(nrow(wc_s[wc_s$stop>time[j],]))}
			prob.pvs <- ifelse(prob.pvs>1,1,prob.pvs)
			prob.pvs <- as.matrix(prob.pvs)
			w[i] <- 1/(1 - exp(t(log(1-prob.pvs))%*%as.matrix(rep(1,nrow(prob.pvs)))))}}
			

#Indicator for duplicating records, 1 if both controls, 2 if 1 case the other control
ncc_o <- ncc[order(ncc$TWINNR,ncc$cc),]

dupl<-NULL
for (i in 2:nrow(ncc_o)-1){
	if (ncc_o$TWINNR[i]==ncc_o$TWINNR[i+1] & ncc_o$cc[i]==ncc_o$cc[i+1]){ncc_o$dupl[i]<-1}
	else if (ncc_o$TWINNR[i]==ncc_o$TWINNR[i+1] & ncc_o$cc[i]!=ncc_o$cc[i+1]){ncc_o$dupl[i]<-2}
	else {ncc_o$dupl[i]<-0}}

#Dataset with duplication indicator for analysis
ncc <- ncc_o[order(ncc_o$rstime,ncc_o$cc,ncc_o$stop),]
	
#Dataset with no duplicated subjects
Uncc <- ncc[ncc$dupl<1,]

# We scale weight so that sum=total population and we eliminate duplicate records
div <- sum(w[ncc$dupl<1 & ncc$cc==0])/nrow(wc_s[wc_s$CVD==0,])
w <- ifelse(w[ncc$dupl<1]==1,1,w[ncc$dupl<1]/div)

#Excluding duplicated records
Rt_ncc_wb_m <- Rt_ncc_wb[ncc$dupl<1]



			
#GOF
ncc$cutRt <- as.factor(cut(log(Rt_ncc_wb), quantile(log(Rt_ncc_wb),(0:5)/5), labels=1:5, include.lowest=T))

res_ncc <- as.numeric(residuals(NCC, type=c("martingale")))
sub_ncc <- ncc$cc-res_ncc

pred_ncc <- aggregate(sub_ncc,list(ncc$cutRt),sum)

obs_ncc <- aggregate(ncc$cc,list(ncc$cutRt),sum)

N_ncc <- table(ncc$cutRt)

HL_ncc <- sum ((obs_ncc$x- pred_ncc$x)^2/(pred_ncc$x))

#GOF
ncc$cutRt2 <- as.factor(cut(Rt_ncc_wb, c(0,0.05,0.2,1), labels=1:3, include.lowest=T))

res_ncc2 <- as.numeric(residuals(NCC, type=c("martingale")))
sub_ncc2 <- ncc$cc-res_ncc2

pred_ncc2 <- aggregate(sub_ncc2,list(ncc$cutRt2),sum)

obs_ncc2 <- aggregate(ncc$cc,list(ncc$cutRt2),sum)

if (sum(as.numeric(pred_ncc2$Group.1))!=6){
	pred_ncc2 <- rbind(pred_ncc2,c(3,0))
	obs_ncc2 <- rbind(obs_ncc2,c(3,0))}



N_ncc2 <- table(ncc$cutRt2)

HL_ncc2 <- sum ((obs_ncc2$x- pred_ncc2$x)^2/(pred_ncc2$x))


ncc$cutRt3 <- as.factor(cut(log(Rt_ncc_wb), wtd.quantile(log(Rt_ncc_wb_m),weights=w,(0:5)/5), labels=1:5, include.lowest=T))
ncc$cutRt3 <- ifelse(log(Rt_ncc_wb)>wtd.quantile(log(Rt_ncc_wb_m),weights=w,(0:5)/5)[6],5,
               ifelse(log(Rt_ncc_wb)<wtd.quantile(log(Rt_ncc_wb_m),weights=w,(0:5)/5)[1],1,ncc$cutRt3))

res_ncc3 <- as.numeric(residuals(NCC, type=c("martingale")))
sub_ncc3 <- ncc$cc-res_ncc3

pred_ncc3 <- aggregate(sub_ncc3,list(ncc$cutRt3),sum)

obs_ncc3 <- aggregate(ncc$cc,list(ncc$cutRt3),sum)

N_ncc3 <- table(ncc$cutRt3)

HL_ncc3 <- sum ((obs_ncc3$x- pred_ncc3$x)^2/(pred_ncc3$x))



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

st0 <-1
st1 <- length(wc_s$bio[wc_s$strata==1 & wc_s$CVD==0])/length(wc_s$bio[wc_s$strata==1 & wc_s$ind==1])
st2 <- length(wc_s$bio[wc_s$strata==2 & wc_s$CVD==0])/length(wc_s$bio[wc_s$strata==2 & wc_s$ind==1])
st3 <- length(wc_s$bio[wc_s$strata==3 & wc_s$CVD==0])/length(wc_s$bio[wc_s$strata==3 & wc_s$ind==1])
st4 <- length(wc_s$bio[wc_s$strata==4 & wc_s$CVD==0])/length(wc_s$bio[wc_s$strata==4 & wc_s$ind==1])

#Redefined weights (cases excluded)
w <- ifelse(cc$CVD==1,st0,
     ifelse(cc$strata==1,st1,
     ifelse(cc$strata==2,st2,
     ifelse(cc$strata==3,st3,st4))))
 
#Excluding cases in the subcohort (The same because we worked on cc)
Rt_BII_wb_m <- Rt_BII_wb

exp <- NULL
for (i in 1:nrow(dstrat2)){
	if (sum(H0_BII$time<=dstrat2$phase1$sample$variables$stop[i])!=0){
	St <- H0_BII$hazard[which(H0_BII$time==max(H0_BII$time[H0_BII$time<=dstrat2$phase1$sample$variables$stop[i]]))] * (1/dstrat2$prob[i])}
	else {St <- 0 }

	temp <- (St*expXB_BII[i])
	exp <- c(exp,temp)
	}	

		

#GOF
dstrat2$phase1$sample$variables$cutY <- as.factor(cut(log(Rt_BII_wb), quantile(log(Rt_BII_wb),(0:5)/5), labels=1:5, include.lowest=T))


res_BII <- dstrat2$phase1$sample$variables$CVD-exp
sub_BII <- exp 

pred_BII <- aggregate(sub_BII,list(dstrat2$phase1$sample$variables$cutY),sum)

obs_BII <- aggregate(dstrat2$phase1$sample$variables$CVD,list(dstrat2$phase1$sample$variables$cutY),sum)

N_BII <- table(dstrat2$phase1$sample$variables$cutY)


HL_BII <- sum ((obs_BII$x-pred_BII$x)^2/(pred_BII$x))


#GOF
dstrat2$phase1$sample$variables$cutY2 <- as.factor(cut(Rt_BII_wb, c(0,0.5,0.20,1), labels=1:3, include.lowest=T))


res_BII2 <- dstrat2$phase1$sample$variables$CVD-exp
sub_BII2 <- exp 

pred_BII2 <- aggregate(sub_BII2,list(dstrat2$phase1$sample$variables$cutY2),sum)

obs_BII2 <- aggregate(dstrat2$phase1$sample$variables$CVD,list(dstrat2$phase1$sample$variables$cutY2),sum)


if (sum(as.numeric(pred_BII2$Group.1))!=6){
	pred_BII2 <- rbind(pred_BII2,c(3,0))
	obs_BII2 <- rbind(obs_BII2,c(3,0))}


N_BII2 <- table(dstrat2$phase1$sample$variables$cutY2)


HL_BII2 <- sum ((obs_BII2$x-pred_BII2$x)^2/(pred_BII2$x))


#GOF
dstrat2$phase1$sample$variables$cutY3 <- as.factor(cut(log(Rt_BII_wb), wtd.quantile(log(Rt_BII_wb_m),weights=w,(0:5)/5), labels=1:5, include.lowest=T))
dstrat2$phase1$sample$variables$cutY3 <- ifelse(log(Rt_BII_wb)>wtd.quantile(log(Rt_BII_wb_m),weights=w,(0:5)/5)[6],5,
               ifelse(log(Rt_BII_wb)<wtd.quantile(log(Rt_BII_wb_m),weights=w,(0:5)/5)[1],1,dstrat2$phase1$sample$variables$cutY3))

res_BII3 <- dstrat2$phase1$sample$variables$CVD-exp
sub_BII3 <- exp 

pred_BII3 <- aggregate(sub_BII3,list(dstrat2$phase1$sample$variables$cutY3),sum)

obs_BII3 <- aggregate(dstrat2$phase1$sample$variables$CVD,list(dstrat2$phase1$sample$variables$cutY3),sum)


N_BII3 <- table(dstrat2$phase1$sample$variables$cutY3)


HL_BII3 <- sum ((obs_BII3$x-pred_BII3$x)^2/(pred_BII3$x))




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


  
#Weights
## We reduce the dataset to speed up calculations
temp<- ncc_s[c("stop","cc","age","SEX")]
w=rep(1,nrow(ncc_s))
for(i in 1:nrow(temp)){
	if(temp$cc[i]==0) {
		time <- temp$stop[temp$cc==1 & temp$stop <= temp$stop[i] & temp$age==temp$age[i] & temp$SEX==temp$SEX[i]]
		prob.pvs <- rep(NA,length(time))
		for(j in 1:length(time)) {
			prob.pvs[j] <- (m-1)/(nrow(wc_s[wc_s$stop>time[j] & wc_s$age==temp$age[i] & wc_s$SEX==temp$SEX[i],]))}
			prob.pvs <- ifelse(prob.pvs>1,1,prob.pvs)
			prob.pvs <- as.matrix(prob.pvs)
			w[i] <- 1/(1 - exp(t(log(1-prob.pvs))%*%as.matrix(rep(1,nrow(prob.pvs)))))}}


# Indicator for duplicating records, 1 if both controls, 2 if 1 case the other control
ncc_so <- ncc_s[order(ncc_s$TWINNR,ncc_s$cc),]

dupl<-NULL
for (i in 2:nrow(ncc_so)-1){
	if (ncc_so$TWINNR[i]==ncc_so$TWINNR[i+1] & ncc_so$cc[i]==ncc_so$cc[i+1]){ncc_so$dupl[i]<-1}
	else if (ncc_so$TWINNR[i]==ncc_so$TWINNR[i+1] & ncc_so$cc[i]!=ncc_so$cc[i+1]){ncc_so$dupl[i]<-2}
	else {ncc_so$dupl[i]<-0}}

# Dataset with duplication indicator for analysis
ncc_s <- ncc_so[order(ncc_so$rstime,ncc_so$cc,ncc_so$stop),]
	
# Dataset with no duplicated subjects
Uncc_s <- ncc_s[ncc_s$dupl<1,]


# We scale weight so that sum=total population and we eliminate duplicate records
div <- sum(w[ncc_s$dupl<1 & ncc_s$cc==0])/nrow(wc_s[wc_s$CVD==0,])
w <- ifelse(w[ncc_s$dupl<1]==1,1,w[ncc_s$dupl<1]/div)

    
#Excluding duplicated records
Rt_ncc_s_wb_m <- Rt_ncc_s_wb[ncc_s$dupl<1]




#GOF
ncc_s$cutRt <- as.factor(cut(log(Rt_ncc_s_wb), quantile(log(Rt_ncc_s_wb),(0:5)/5), labels=1:5, include.lowest=T))


res_ncc_s <- as.numeric(residuals(NCC_s, type=c("martingale")))
sub_ncc_s <- ncc_s$cc-res_ncc_s


pred_ncc_s <- aggregate(sub_ncc_s,list(ncc_s$cutRt),sum)

obs_ncc_s <- aggregate(ncc_s$cc,list(ncc_s$cutRt),sum)

N_ncc_s <- table(ncc_s$cutRt)

HL_ncc_s <- sum ((obs_ncc_s$x- pred_ncc_s$x)^2/(pred_ncc_s$x))
 
 
#GOF
ncc_s$cutRt2 <- as.factor(cut(Rt_ncc_s_wb, c(0,0.05,0.2,1), labels=1:3, include.lowest=T))



res_ncc_s2 <- as.numeric(residuals(NCC_s, type=c("martingale")))
sub_ncc_s2 <- ncc_s$cc-res_ncc_s2


pred_ncc_s2 <- aggregate(sub_ncc_s2,list(ncc_s$cutRt2),sum)

obs_ncc_s2 <- aggregate(ncc_s$cc,list(ncc_s$cutRt2),sum)



if (sum(as.numeric(pred_ncc_s2$Group.1))!=6){
	pred_ncc_s2 <- rbind(pred_ncc_s2,c(3,0))
	obs_ncc_s2 <- rbind(obs_ncc_s2,c(3,0))}



N_ncc_s2 <- table(ncc_s$cutRt2)

HL_ncc_s2 <- sum ((obs_ncc_s2$x- pred_ncc_s2$x)^2/(pred_ncc_s2$x))
 

ncc_s$cutRt3 <- as.factor(cut(log(Rt_ncc_s_wb), wtd.quantile(log(Rt_ncc_s_wb_m),weights=w,(0:5)/5), labels=1:5, include.lowest=T))
ncc_s$cutRt3 <- ifelse(log(Rt_ncc_s_wb)>wtd.quantile(log(Rt_ncc_s_wb_m),weights=w,(0:5)/5)[6],5,
               ifelse(log(Rt_ncc_s_wb)<wtd.quantile(log(Rt_ncc_s_wb_m),weights=w,(0:5)/5)[1],1,ncc_s$cutRt3))

res_ncc_s3 <- as.numeric(residuals(NCC_s, type=c("martingale")))
sub_ncc_s3 <- ncc_s$cc-res_ncc_s3

pred_ncc_s3 <- aggregate(sub_ncc_s3,list(ncc_s$cutRt3),sum)

obs_ncc_s3 <- aggregate(ncc_s$cc,list(ncc_s$cutRt3),sum)

N_ncc_s3 <- table(ncc_s$cutRt3)

HL_ncc_s3 <- sum ((obs_ncc_s3$x- pred_ncc_s3$x)^2/(pred_ncc_s3$x))




     
#################
#MERGING RESULTS#
#################

### a.ASSOCIATION ###




### c.PREDICTION MEASURES ###

pred <- cbind(pred_wc$x,pred_Pr$x,pred_ncc$x,pred_BII$x,pred_ncc_s$x,c(1,2,3,4,5))
obs <- cbind(obs_wc$x,obs_Pr$x,obs_ncc$x,obs_BII$x,obs_ncc_s$x,c(1,2,3,4,5))
HL <- cbind(HL_wc,HL_wc_b,HL_Pr,HL_ncc,HL_BII,HL_ncc_s)

pred2 <- cbind(pred_wc2$x,pred_Pr2$x,pred_ncc2$x,pred_BII2$x,pred_ncc_s2$x,c(1,2,3))
obs2 <- cbind(obs_wc2$x,obs_Pr2$x,obs_ncc2$x,obs_BII2$x,obs_ncc_s2$x,c(1,2,3))
HL2 <- cbind(HL_wc2,HL_wc2_b,HL_Pr2,HL_ncc2,HL_BII2,HL_ncc_s2)

pred3 <- cbind(pred_wc$x,pred_Pr3$x,pred_ncc3$x,pred_BII3$x,pred_ncc_s3$x,c(1,2,3,4,5))
obs3 <- cbind(obs_wc$x,obs_Pr3$x,obs_ncc3$x,obs_BII3$x,obs_ncc_s3$x,c(1,2,3,4,5))
HL3 <- cbind(HL_wc,HL_wc_b,HL_Pr3,HL_ncc3,HL_BII3,HL_ncc_s3)

return(list(pred,obs,HL,pred2,obs2,HL2,pred3,obs3,HL3))


}


F11 <-NULL
F11_GOF <- NULL
F11_GOF2 <- NULL
F11_pred <- NULL
F11_pred2 <- NULL
F11_pred3 <- NULL
F11_obs <- NULL
F11_obs2 <- NULL
F11_obs3 <- NULL
F11_HL <- NULL
F11_HL2 <- NULL
F11_HL3 <- NULL



F13 <-NULL
F13_GOF <- NULL
F13_GOF2 <- NULL
F13_pred <- NULL
F13_pred2 <- NULL
F13_pred3 <- NULL
F13_obs <- NULL
F13_obs2 <- NULL
F13_obs3 <- NULL
F13_HL <- NULL
F13_HL2 <- NULL 
F13_HL3 <- NULL 

for (k in 1:nrep) {   F11<- ccfun(229,2,1095)
	                  F11_pred <- rbind(F11_pred, F11[[1]])
	                  F11_obs <- rbind(F11_obs, F11[[2]])
	                  F11_HL <- rbind(F11_HL, F11[[3]])
	                  F11_pred2 <- rbind(F11_pred2, F11[[4]])
	                  F11_obs2 <- rbind(F11_obs2, F11[[5]])
	                  F11_HL2 <- rbind(F11_HL2, F11[[6]])
	                  F11_pred3 <- rbind(F11_pred3, F11[[7]])
	                  F11_obs3 <- rbind(F11_obs3, F11[[8]])
	                  F11_HL3 <- rbind(F11_HL3, F11[[9]])
                      print(k)}



save <- list(F11_pred,F11_obs,F11_HL,F11_pred2,F11_obs2,F11_HL2,F11_pred3,F11_obs3,F11_HL3)
save(save,file="GOF_bio1_1_HL.Rdata")


for (k in 1:nrep) {   F13<- ccfun(632,4,1095)
	                  F13_pred <- rbind(F13_pred, F13[[1]])
	                  F13_obs <- rbind(F13_obs, F13[[2]])
	                  F13_HL <- rbind(F13_HL, F13[[3]])
	                  F13_pred2 <- rbind(F13_pred2, F13[[4]])
	                  F13_obs2 <- rbind(F13_obs2, F13[[5]])
	                  F13_HL2 <- rbind(F13_HL2, F13[[6]])
	                  F13_pred3 <- rbind(F13_pred3, F13[[7]])
	                  F13_obs3 <- rbind(F13_obs3, F13[[8]])
	                  F13_HL3 <- rbind(F13_HL3, F13[[9]])
	                                      print(k)}



save <- list(F13_pred,F13_obs,F13_HL,F13_pred2,F13_obs2,F13_HL2,F13_pred3,F13_obs3,F13_HL3)
save(save,file="GOF_bio1_3_HL.Rdata")



	