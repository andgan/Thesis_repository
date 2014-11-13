#----------------------------------------------
# Filename: All_paper.R
# Study: CCNCC
# Author: Andrea Ganna
# Date: 08OCT2010
# Updated: 15DEC2010
# Purpose: All main analyses for CCNCC
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



nrep <- 200


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


### LOAD FUNCTIONS ###

# HOSMER-LAMESHOW
HL=function(data,pred,out,w){
des <-svydesign(id=~1, weights=w, data=data)
breaks <- svyquantile(pred, des, (0:10)/10)
cutRt <-cut(pred, breaks, labels=1:(length(breaks)-1), include.lowest=T)
obs <- xtabs (cbind((1-out)*w, out) ~ cutRt)
wm<-NULL
wmm<-NULL
for (i in 1:(length(breaks)-1)){
wm<-weighted.mean(pred[cutRt==i],w[cutRt==i])
wmm <- c(wmm,wm)}
tot<-as.numeric(obs[,1]+obs[,2])
expected_e<-wmm*tot
expected<-cbind(tot-expected_e,expected_e)
HL <- sum ((obs[,2]- expected[,2])^2/(expected[,2]*(1-wmm)))
pHL <- 1 - pchisq (HL, (length(breaks)-1)-2)
return(data.frame(HL))}


# HOSMER-LAMESHOW FAST
HL_fast=function(pred,out){
cutRt <-cut(pred, quantile(pred,(0:10)/10), labels=1:10, include.lowest=T)
obs <- as.numeric(xtabs ( out ~ cutRt))
wmm <- aggregate(pred,by=list(cutRt),mean)
tot<-as.numeric(table(cutRt))
expected<-wmm$x*tot
HL <- sum ((obs- expected)^2/(expected*(1-wmm$x)))
return(data.frame(HL))}



#NRI
NRI=function(wc_,wcout,predwith,predwithout,out,cat,w){
P1 <- predwith
P2 <- predwithout
cutP1_e <-cut(P1[out==1], breaks=cat, label=1:3)
cutP1_ne <-cut(P1[out==0], breaks=cat, label=1:3)
cutP2_e <-cut(P2[out==1], breaks=cat, label=1:3)
cutP2_ne <-cut(P2[out==0], breaks=cat, label=1:3)
# NON EVENTS
down_ne <- NULL
up_ne <- NULL
dow<-NULL
up <- NULL
for (i in 1:length(cutP1_ne)){
	if (as.numeric(cutP2_ne[i]) > as.numeric(cutP1_ne[i])){
		dow = w[i]}
	else {dow=0}
	if  (as.numeric(cutP2_ne[i]) < as.numeric(cutP1_ne[i])){
		up = w[i]}
	else {up=0}
down_ne <- sum(c(down_ne,dow))
up_ne <- sum(c(up_ne,up))}
down_ne <-down_ne/nrow(wc_[wcout==0,])
up_ne <- up_ne/nrow(wc_[wcout==0,])
#EVENTS
down_e <- NULL
up_e <- NULL
dow<-NULL
up <- NULL
for (i in 1:length(cutP1_e)){
	if (as.numeric(cutP2_e[i]) > as.numeric(cutP1_e[i])){
		dow = 1}
	else {dow=0}
	if  (as.numeric(cutP2_e[i]) < as.numeric(cutP1_e[i])){
		up = 1}
	else {up=0}
down_e <- sum(c(down_e,dow))
up_e <- sum(c(up_e,up))}
down_e <-down_e/nrow(wc_[wcout==1,])
up_e <- up_e/nrow(wc_[wcout==1,])
NRI <- ((up_e-down_e)-(up_ne-down_ne))*100
seNRI <- ((up_e+down_e)/length(cutP1_e)+(up_ne+down_ne)/length(cutP1_ne))^0.5
zNRI <- (NRI/100)/seNRI
pNRI <- round(2 * (pnorm(-abs(zNRI))),2)
return(data.frame(NRI,seNRI))}


#IDI
IDI=function(predwith,predwithout,out,w){
P1 <- predwith
P2 <- predwithout
IDI <- (weighted.mean(P1[out==1],w[out==1])-weighted.mean(P2[out==1],w[out==1]))-(weighted.mean(P1[out==0],w[out==0])-weighted.mean(P2[out==0],w[out==0]))
d <- P1-P2
seIDI <- ((sd(d[out==0])/length(d[out==0])^0.5)^2+(sd(d[out==1])/length(d[out==1])^0.5)^2)^0.5
zIDI <- IDI/seIDI
pIDI <- round(2 * (pnorm(-abs(zIDI))),2)
return(data.frame(IDI,seIDI))}


#C-INDEX
CIND=function(pred,out,time,id,w){
    sumci<- NULL
	sumdi<- NULL
	ci <- NULL
	di <- NULL
for (i in 1:length(id)){
	if (out[i]==1){
		ci <- sum(w[((time[i]>time & pred>pred[i]) | (time>time[i] & pred[i]>pred)) & out==1 & id[i]!=id]) + 2*sum(w[(time>=time[i] & pred[i]>pred) & out==0 & id[i]!=id])
	    di <- sum(w[((time[i]>time & pred<pred[i]) | (time>time[i] & pred[i]<pred)) & out==1 & id[i]!=id]) + 2*sum(w[(time>=time[i] & pred[i]<pred) & out==0 & id[i]!=id])}
	else{
		ci<-0
		di<-0}
	sumci <- sum(c(sumci,ci))
	sumdi <- sum(c(sumdi,di))}
	cindex <- sumci/(sumci+sumdi)
	return(data.frame(cindex))}



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



#####


#Without Biomarker
GOLD <-coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp  + diabetess, data=wc_all)

#XB
cf <- GOLD$coefficient
expXB_go <- exp(cf[1]*wc_all$age + cf[2]*wc_all$TC + cf[3]*wc_all$SEX  + cf[4]*wc_all$smoke + cf[5]*wc_all$antihyp + cf[6]*wc_all$sbp + cf[7]*+wc_all$diabetess)

#Baseline Hazard
H0_go <- basehaz(GOLD, center=FALSE)

#Find the t-years H0	
St_go <- exp(-H0_go$hazard[which(abs(H0_go$time-t)==min(abs(H0_go$time-t)))])
	
#For prediction measure
Rt_go_wob <- (1-St_go^expXB_go)



### c.PREDICTION MEASURES ###


#GOF
wc_all$cutRt <- as.factor(cut(log(Rt_go_wb), quantile(log(Rt_go_wb),(0:5)/5), labels=1:5, include.lowest=T))
GOLD2 <- coxph(Surv(start,stop,CVD)~age + TC  + SEX  +  smoke + antihyp +sbp + bio + diabetess + cutRt, data=wc_all)

GOF<- try(wald.test(vcov(GOLD2),coef(GOLD2), Terms=9:12),silent=T)
GOF_go_wb <- try(cbind(as.vector(GOF$result$chi2)[1],"GOF_go_wb"),silent=T)
colnames(GOF_go_wb) <- c("GOF","name")

#H-L FAST
HL_go_wb <- cbind(HL_fast(Rt_go_wb,wc_all$CVD),"HL_go_wb")
colnames(HL_go_wb) <- c("HL","name")


#NRI        
NRI_go <- cbind(NRI(wc_all, wc_all$CVD, Rt_go_wb,Rt_go_wob,wc_all$CVD, c(0,0.05,0.20,100), matrix(1,ncol=1, nrow=nrow(wc_all[wc_all$CVD==0,]))),"NRI_go")
colnames(NRI_go) <- c("NRI","name")

#IDI
IDI_go <- cbind(IDI(Rt_go_wb,Rt_go_wob,wc_all$CVD,matrix(1,ncol=1, nrow=nrow(wc_all))),"IDI_go")
colnames(IDI_go) <- c("IDI","name")

#CINDEX	
CIND_go_wb<-cbind(CIND(Rt_go_wb,wc_all$CVD,wc_all$stop,wc_all$TWINNR, matrix(1,ncol=1, nrow=nrow(wc_all))),"CIND_go_wb")
CIND_go_wob<-cbind(CIND(Rt_go_wob,wc_all$CVD,wc_all$stop,wc_all$TWINNR, matrix(1,ncol=1, nrow=nrow(wc_all))),"CIND_go_wob")
colnames(CIND_go_wob) <- c("cindex","name")
colnames(CIND_go_wb) <- c("cindex","name")




 
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


#Find the t-years H0	
St_wc <- exp(-H0_wc$hazard[which(abs(H0_wc$time-t)==min(abs(H0_wc$time-t)))])
	
#Individual risk
Rt_wc <- (1-St_wc^expXB_wc)*100

Rt_wc <- data.frame(cbind(Rt_wc,wc_s$TWINNR_or))
colnames(Rt_wc) <- c("Risk at t-years","TWINNR")

#For prediction measure
Rt_wc_wb <- (1-St_wc^expXB_wc)




#####


#Without Biomarker
WC <-coxph(Surv(start,stop,CVD)~age + TC  + SEX + smoke + antihyp +sbp + diabetess, data=wc_s)


#XB
cf <- WC$coefficient
expXB_wc <- exp(cf[1]*wc_s$age + cf[2]*wc_s$TC + cf[3]*wc_s$SEX +  cf[4]*wc_s$smoke + cf[5]*wc_s$antihyp + cf[6]*wc_s$sbp + cf[7]*+wc_s$diabetess)

#Baseline Hazard
H0_wc <- basehaz(WC, center=FALSE)
H0_wc2 <- basehaz(WC)

#Find the t-years H0
St_wc <- exp(-H0_wc$hazard[which(abs(H0_wc$time-t)==min(abs(H0_wc$time-t)))])
St_wc_wob <- exp(-H0_wc2$hazard[which(abs(H0_wc2$time-t)==min(abs(H0_wc2$time-t)))])
	
#For prediction measure
Rt_wc_wob <- (1-St_wc^expXB_wc)



#Without Biomarker,SEX and age
WC <-coxph(Surv(start,stop,CVD)~ TC  + smoke + antihyp +sbp + diabetess, data=wc_s)

#Baseline Hazard
H0_wc <- basehaz(WC)

#Find the t-years H0
St_wc_wob3 <- exp(-H0_wc$hazard[which(abs(H0_wc$time-t)==min(abs(H0_wc$time-t)))])



### c.PREDICTION MEASURES ###

#We break individual risk at random
dupl <-  duplicated(Rt_wc_wob)
random<-rnorm(length(Rt_wc_wob),0,1)/10000
Rt_wc_wob <- ifelse(dupl,Rt_wc_wob+random,Rt_wc_wob)

dupl <-  duplicated(Rt_wc_wb)
random<-rnorm(length(Rt_wc_wb),0,1)/10000
Rt_wc_wb <- ifelse(dupl,Rt_wc_wb+random,Rt_wc_wb)



#GOF
wc_s$cutRt <- as.factor(cut(log(Rt_wc_wb), quantile(log(Rt_wc_wb),(0:5)/5), labels=1:5, include.lowest=T))
tag <- 0
withCallingHandlers(
WC2 <- coxph(Surv(start,stop,CVD)~age + TC  + SEX  + smoke + antihyp +sbp + bio + diabetess + cutRt, data=wc_s),
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(u <- vcov(WC2),coef(WC2), Terms=9:12)
GOF_wc_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_wc_wb")
colnames(GOF_wc_wb) <- c("GOF","name")} else {GOF_wc_wb <- c(NA,"GOF_wc_wb")}




#H-L FAST
HL_wc_wb <- cbind(HL_fast(Rt_wc_wb,wc_s$CVD),"HL_wc_wb")
colnames(HL_wc_wb) <- c("HL","name")


#NRI        
NRI_wc <- cbind(NRI(wc_s, wc_s$CVD, Rt_wc_wb,Rt_wc_wob,wc_s$CVD, c(0,0.05,0.20,100), matrix(1,ncol=1, nrow=nrow(wc_s[wc_s$CVD==0,]))),"NRI_wc")
colnames(NRI_wc) <- c("NRI","name")

#IDI
IDI_wc <- cbind(IDI(Rt_wc_wb,Rt_wc_wob,wc_s$CVD,matrix(1,ncol=1, nrow=nrow(wc_s))),"IDI_wc")
colnames(IDI_wc) <- c("IDI","name")

#CINDEX	
CIND_wc_wb<-cbind(CIND(Rt_wc_wb,wc_s$CVD,wc_s$stop,wc_s$TWINNR, matrix(1,ncol=1, nrow=nrow(wc_s))),"CIND_wc_wb")
CIND_wc_wob<-cbind(CIND(Rt_wc_wob,wc_s$CVD,wc_s$stop,wc_s$TWINNR, matrix(1,ncol=1, nrow=nrow(wc_s))),"CIND_wc_wob")
colnames(CIND_wc_wob) <- c("cindex","name")
colnames(CIND_wc_wb) <- c("cindex","name")




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

Rt_Pr <- data.frame(cbind(Rt_Pr,cc_s$TWINNR_or,cc_s$CVD_cc))
colnames(Rt_Pr) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_Pr_wb <- (1-St_Pr^cc_s$expXB_Pr)

#Risk Excluding cases in the subcohort
Rt_Pr_wb_m <- Rt_Pr_wb[(cc_s$CVD_cc==0 & cc_s$CVD==1)==F]



#####

#Individual risk with H0 from the whole cohort

#XB

cf <- PR$coefficient
expXB_Pr <- exp(cf[1]*(cc$age-mean(wc$age)) + cf[2]*(cc$TC-mean(wc$TC)) + cf[3]*(cc$SEX-mean(wc$SEX)) + cf[4]*(cc$smoke-mean(wc$smoke)) + cf[5]*(cc$antihyp-mean(wc$antihyp)) + cf[6]*(cc$sbp-mean(wc$sbp)) + cf[7]*(cc$bio-mean(cc$bio[cc$CVD==0])) + cf[8]*(cc$diabetess-mean(wc$diabetess)))

#Individual Risk
Rt_Pr_wc <- (1-St_wc_wob^expXB_Pr)*100

Rt_Pr_wc <- data.frame(cbind(Rt_Pr_wc,cc$TWINNR_or,cc$CVD))
colnames(Rt_Pr_wc) <- c("Risk at t-years","TWINNR")

#####



#Without biomarker

PR <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + diabetess + cluster(as.factor(keys)), data=cc_s)

#XB
cf <- PR$coefficient
cc_s$expXB_Pr <- exp(cf[1]*cc_s$age + cf[2]*cc_s$TC + cf[3]*cc_s$SEX + cf[4]*cc_s$smoke + cf[5]*cc_s$antihyp + cf[6]*cc_s$sbp + cf[7]*cc_s$diabetess )

#Creating failures time (all CVD events)
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

#For prediction measures
Rt_Pr_wob <- (1-St_Pr^cc_s$expXB_Pr)

#Risk Excluding cases in the subcohort
Rt_Pr_wob_m <- Rt_Pr_wob[(cc_s$CVD_cc==0 & cc_s$CVD==1)==F]


### c.PREDICTION MEASURES ###

#Weights for prediction measures
w <- 1/(ifelse(cc$CVD==1,1,nrow(cc[wc_s$ind==1 & wc_s$CVD==0,])/nrow(wc_s[wc_s$CVD==0,])))
 
#We break individual risk at random
dupl <-  duplicated(Rt_Pr_wob_m)
random<-rnorm(length(Rt_Pr_wob_m),0,1)/10000
Rt_Pr_wob_m <- ifelse(dupl,Rt_Pr_wob_m+random,Rt_Pr_wob_m)

dupl <-  duplicated(Rt_Pr_wb_m)
random<-rnorm(length(Rt_Pr_wb_m),0,1)/10000
Rt_Pr_wb_m <- ifelse(dupl,Rt_Pr_wb_m+random,Rt_Pr_wb_m)


#GOF
cc_s$cutRt <- as.factor(cut(log(Rt_Pr_wb), quantile(log(Rt_Pr_wb),(0:5)/5), labels=1:5, include.lowest=T))
tag <- 0
withCallingHandlers(
PR2 <- coxph(Surv(start_cc,stop_cc,CVD_cc)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess +   cutRt + cluster(as.factor(keys)), data=cc_s),
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(vcov(PR2),coef(PR2), Terms=9:12)
GOF_Pr_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_Pr_wb")
colnames(GOF_Pr_wb) <- c("GOF","name")} else {GOF_Pr_wb <- c(NA,"GOF_Pr_wb")}


#H-L
HL_Pr_wb <- cbind(HL(cc,Rt_Pr_wb_m,cc$CVD,w),"HL_Pr_wb")
colnames(HL_Pr_wb) <- c("HL","name")

#NRI 
NRI_Pr <- cbind(NRI(wc_s, wc_s$CVD, Rt_Pr_wb_m,Rt_Pr_wob_m,cc$CVD, c(0,0.05,0.20,100), w[cc$CVD==0]),"NRI_Pr")
colnames(NRI_Pr) <- c("NRI","name")


#IDI 
IDI_Pr <- cbind(IDI(Rt_Pr_wb_m,Rt_Pr_wob_m,cc$CVD,w),"IDI_Pr")
colnames(IDI_Pr) <- c("IDI","name")


#C-INDEX 
CIND_Pr_wb<-cbind(CIND(Rt_Pr_wb_m,cc$CVD,cc$stop,cc$TWINNR,w),"CIND_Pr_wb")
colnames(CIND_Pr_wb) <- c("cindex","name")
CIND_Pr_wob<-cbind(CIND(Rt_Pr_wob_m,cc$CVD,cc$stop,cc$TWINNR,w),"CIND_Pr_wob")
colnames(CIND_Pr_wob) <- c("cindex","name")



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
Rt_ncc <- data.frame(cbind(Rt_ncc,ncc$TWINNR_or,ncc$cc))
colnames(Rt_ncc) <- c("Risk at t-years","TWINNR")

# For prediction measures
Rt_ncc_wb <- (1-St_ncc^expXB_ncc)


#####

#Individual risk with H0 from the wholw cohort

#XB

#XB 
cf <- NCC$coefficient
expXB_ncc <- exp(cf[1]*(ncc$age-mean(wc$age)) + cf[2]*(ncc$TC-mean(wc$TC)) + cf[3]*(ncc$SEX-mean(wc$SEX)) + cf[4]*(ncc$smoke-mean(wc$smoke)) + cf[5]*(ncc$antihyp-mean(wc$antihyp)) + cf[6]*(ncc$sbp-mean(wc$sbp)) + cf[7]*(ncc$bio-mean(ncc$bio[ncc$cc==0])) + cf[8]*(ncc$diabetess-mean(wc$diabetess)))

#Individual risk	
Rt_ncc_wc <- (1-St_wc_wob^expXB_ncc)*100
Rt_ncc_wc <- data.frame(cbind(Rt_ncc_wc,ncc$TWINNR_or,ncc$cc))
colnames(Rt_ncc_wc) <- c("Risk at t-years","TWINNR")




#####


#Without biomarker 
NCC <- coxph(Surv(fakentry,rstime,cc)~age+TC+ SEX + smoke + antihyp + sbp + diabetess + offset(logw), data=ncc) 

#XB  
cf <- NCC$coefficient
expXB_ncc <- exp(cf[1]*ncc$age + cf[2]*ncc$TC + cf[3]*ncc$SEX  + cf[4]*ncc$smoke + cf[5]*ncc$antihyp + cf[6]*ncc$sbp  +cf[7]*ncc$diabetess)

#Baseline Hazard
H0_ncc <- basehaz(NCC, center=FALSE)

#Find the t-years H0wc	
St_ncc <- exp(-H0_ncc$hazard[which(abs(H0_ncc$time-t)==min(abs(H0_ncc$time-t)))])

# For prediction measures
Rt_ncc_wob <- (1-St_ncc^expXB_ncc)


### c.PREDICTION MEASURES ###


#Weights and duplicated data for prediction measures

#Weights

# We reduce the dataset to speed up calculations
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

#Excluding duplicated records
Rt_ncc_wob_m <- Rt_ncc_wob[ncc$dupl<1]


# Old Weights
#w<-ifelse(Uncc$CVD==1,1,nrow(wc_s[wc_s$CVD==0,])/nrow(Uncc[Uncc$CVD==0,]))

#We break individual risk at random
dupl <-  duplicated(Rt_ncc_wob_m)
random<-rnorm(length(Rt_ncc_wob_m),0,1)/10000
Rt_ncc_wob_m <- ifelse(dupl,Rt_ncc_wob_m+random,Rt_ncc_wob_m)

dupl <-  duplicated(Rt_ncc_wb_m)
random<-rnorm(length(Rt_ncc_wb_m),0,1)/10000
Rt_ncc_wb_m  <- ifelse(dupl,Rt_ncc_wb_m +random,Rt_ncc_wb_m )

#GOF
ncc$cutRt <- as.factor(cut(log(Rt_ncc_wb), quantile(log(Rt_ncc_wb),(0:5)/5), labels=1:5, include.lowest=T))
tag <- 0
withCallingHandlers(
NCC2 <- coxph(Surv(fakentry,rstime,cc)~age+TC+ SEX + smoke + antihyp + sbp + bio + diabetess + cutRt + offset(logw), data=ncc), 
warning=function(w) {tag <<-1})

if(tag==0) {GOF<- wald.test(vcov(NCC2),coef(NCC2), Terms=9:12)
GOF_ncc_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_ncc_wb")
colnames(GOF_ncc_wb) <- c("GOF","name")} else {GOF_ncc_wb <- c(NA,"GOF_ncc_wb")}



#H-L
HL_ncc_wb <- cbind(HL(Uncc,Rt_ncc_wb_m,Uncc$cc,w),"HL_ncc_wb")
colnames(HL_ncc_wb) <- c("HL","name")

#NRI 
NRI_ncc <- cbind(NRI(wc_s, wc_s$CVD, Rt_ncc_wb_m,Rt_ncc_wob_m,Uncc$cc, c(0,0.05,0.20,100), w[Uncc$cc==0]),"NRI_ncc")
colnames(NRI_ncc) <- c("NRI","name")

#IDI 
IDI_ncc <- cbind(IDI(Rt_ncc_wb_m,Rt_ncc_wob_m,Uncc$cc,w),"IDI_ncc")
colnames(IDI_ncc) <- c("IDI","name")


#C-INDEX 
CIND_ncc_wb<-cbind(CIND(Rt_ncc_wb_m,Uncc$cc,Uncc$stop,Uncc$TWINNR,w),"CIND_ncc_wb")
CIND_ncc_wob<-cbind(CIND(Rt_ncc_wob_m,Uncc$cc,Uncc$stop,Uncc$TWINNR,w),"CIND_ncc_wob")
colnames(CIND_ncc_wob) <- c("cindex","name")
colnames(CIND_ncc_wb) <- c("cindex","name")




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
Rt_BII <- data.frame(cbind(Rt_BII,cc$TWINNR_or,cc$CVD))
colnames(Rt_BII) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_BII_wb <- (1-St_BII^expXB_BII)

#####

#Individual risk with H0 from the wholw cohort

#XB

cf <- BII_s4$coefficient
expXB_BII <- exp(cf[1]*(cc$age-mean(wc$age)) + cf[2]*(cc$TC-mean(wc$TC)) + cf[3]*(cc$SEX-mean(wc$SEX))  + cf[4]*(cc$smoke-mean(wc$smoke)) + cf[5]*(cc$antihyp-mean(wc$antihyp)) + cf[6]*(cc$sbp-mean(wc$sbp)) + cf[7]*(cc$bio-mean(cc$bio[cc$CVD==0])) + cf[8]*(cc$diabetess-mean(wc$diabetess)))

	
#Individual risk	
Rt_BII_wc <- (1-St_wc_wob^expXB_BII)*100
Rt_BII_wc <- data.frame(cbind(Rt_BII_wc,cc$TWINNR_or,cc$CVD))
colnames(Rt_BII_wc) <- c("Risk at t-years","TWINNR")

#####



# Without biomarker
BII_s4<-svycoxph(Surv(start,stop,CVD)~age+TC+ SEX + smoke + antihyp + sbp + diabetess,design=dstrat2)  

#XB
cf <- BII_s4$coefficient
expXB_BII <- exp(cf[1]*cc$age + cf[2]*cc$TC + cf[3]*cc$SEX + cf[4]*cc$smoke + cf[5]*cc$antihyp + cf[6]*cc$sbp + cf[7]*cc$diabetess)

#Baseline Hazard
H0_BII <- basehaz(BII_s4, center=FALSE)

#Find the t-years H0wc	
St_BII <- exp(-H0_BII$hazard[which(abs(H0_BII$time-t)==min(abs(H0_BII$time-t)))])
	
#For prediction
Rt_BII_wob <- (1-St_BII^expXB_BII)


### c.PREDICTION MEASURES ###

#Weights and duplicated data for prediction measures

# Redifing weights including only controls, cases are weighted with w=1
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
Rt_BII_wob_m <- Rt_BII_wob

#We break individual risk at random
dupl <-  duplicated(Rt_BII_wob_m)
random<-rnorm(length(Rt_BII_wob_m),0,1)/10000
Rt_BII_wob_m <- ifelse(dupl,Rt_BII_wob_m+random,Rt_BII_wob_m)

dupl <-  duplicated(Rt_BII_wb_m)
random<-rnorm(length(Rt_BII_wb_m),0,1)/10000
Rt_BII_wb_m  <- ifelse(dupl,Rt_BII_wb_m +random,Rt_BII_wb_m )


#GOF
dstrat2$cutRt <- as.factor(cut(log(Rt_BII_wb), quantile(log(Rt_BII_wb),(0:5)/5), labels=1:5, include.lowest=T))
tag <- 0
withCallingHandlers(
BII_s42 <-svycoxph(Surv(start,stop,CVD)~age+TC+ SEX  + smoke + antihyp + sbp + bio + diabetess + cutRt,design=dstrat2), 
warning=function(w) {tag <<-1}, error=function(e) {tag<<-1})

if(tag==0) {GOF <- wald.test(vcov(BII_s42),coef(BII_s42), Terms=9:12)
GOF_BII_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_BII_wb")
colnames(GOF_BII_wb) <- c("GOF","name")} else {GOF_BII_wb<- c(NA,"GOF_BII_wb")}


        
#H-L
HL_BII_wb <- cbind(HL(cc,Rt_BII_wb_m,cc$CVD,w),"HL_BII_wb")
colnames(HL_BII_wb) <- c("HL","name")


#NRI 
NRI_BII <- cbind(NRI(wc_s, wc_s$CVD, Rt_BII_wb_m,Rt_BII_wob_m,cc$CVD, c(0,0.05,0.20,100), w[cc$CVD==0]),"NRI_BII")
colnames(NRI_BII) <- c("NRI","name")


#IDI 
IDI_BII <- cbind(IDI(Rt_BII_wb_m,Rt_BII_wob_m,cc$CVD,w), "IDI_BII")
colnames(IDI_BII) <- c("IDI","name")


#C-INDEX 
CIND_BII_wb<-cbind(CIND(Rt_BII_wb_m,cc$CVD,cc$stop,cc$TWINNR,w),"CIND_BII_wb")
CIND_BII_wob<-cbind(CIND(Rt_BII_wob_m,cc$CVD,cc$stop,cc$TWINNR,w),"CIND_BII_wob")
colnames(CIND_BII_wb) <- c("cindex","name")
colnames(CIND_BII_wob) <- c("cindex","name")



 
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
Rt_ncc_s <- data.frame(cbind((1-(St_ncc_s^expXB_ncc_s))*100,ncc_s$TWINNR_or,ncc_s$cc))
colnames(Rt_ncc_s) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_ncc_s_wb <- Rt_ncc_s[,1]/100


#####

#Individual risk with H0 from the whole cohort

#XB
cf <- NCC_s$coefficient
expXB_ncc_s <- exp(cf[1]*(ncc_s$TC-mean(wc$TC))  + cf[2]*(ncc_s$smoke-mean(wc$smoke)) + cf[3]*(ncc_s$antihyp-mean(wc$antihyp)) + cf[4]*(ncc_s$sbp-mean(wc$sbp)) + cf[5]*(ncc_s$bio-mean(ncc_s$bio[ncc_s$cc==0])) + cf[6]*(ncc_s$diabetess-mean(wc$diabetess)))



#Individual risk
Rt_ncc_s_wc <- data.frame(cbind((1-(St_wc_wob3^expXB_ncc_s))*100,ncc_s$TWINNR_or,ncc_s$cc))
colnames(Rt_ncc_s_wc) <- c("Risk at t-years","TWINNR")



#####

#Without biomarkers
           
NCC_s <- coxph(Surv(fakentry,rstime,cc)~TC  + smoke + antihyp + sbp + diabetess, data=ncc_s) 

#XB  
cf <- NCC_s$coefficient
expXB_ncc_s <- exp(cf[1]*ncc_s$TC  + cf[2]*ncc_s$smoke + cf[3]*ncc_s$antihyp + cf[4]*ncc_s$sbp +  cf[5]*ncc_s$diabetess)
sexpXB_ncc_s <- aggregate(expXB_ncc_s*(ncc_s$w),list(time=ncc_s$rstime),sum)

#Baseline Hazard
H0_ncc_s<- data.frame(cbind(sexpXB_ncc_s$time,(cumsum(1/sexpXB_ncc_s$x))))
colnames(H0_ncc_s) <- c("time","hazard")

#Find the t-years 
St_ncc_s <- exp(-H0_ncc_s$hazard[which(abs(H0_ncc_s$time-t)==min(abs(H0_ncc_s$time-t)))])

#Individual risk
Rt_ncc_s <- data.frame(cbind((1-(St_ncc_s^expXB_ncc_s))*100,ncc_s$TWINNR_or,ncc_s$cc))
colnames(Rt_ncc_s) <- c("Risk at t-years","TWINNR")

#For prediction measures
Rt_ncc_s_wob <- Rt_ncc_s[,1]/100




### c.PREDICTION MEASURES ###

  
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

#Excluding duplicated records
Rt_ncc_s_wob_m <- Rt_ncc_s_wob[ncc_s$dupl<1]


#We break individual risk at random
dupl <-  duplicated(Rt_ncc_s_wob_m)
random<-rnorm(length(Rt_ncc_s_wob_m),0,1)/10000
Rt_ncc_s_wob_m <- ifelse(dupl,Rt_ncc_s_wob_m+random,Rt_ncc_s_wob_m)

dupl <-  duplicated(Rt_ncc_s_wb_m)
random<-rnorm(length(Rt_ncc_s_wb_m),0,1)/10000
Rt_ncc_s_wb_m <- ifelse(dupl,Rt_ncc_s_wb_m+random,Rt_ncc_s_wb_m)



#GOF
ncc_s$cutRt <- as.factor(cut(log(Rt_ncc_s_wb), quantile(log(Rt_ncc_s_wb),(0:5)/5), labels=1:5, include.lowest=T))
tag <- 0
withCallingHandlers(
NCC_s2 <- coxph(Surv(fakentry,rstime,cc)~TC + smoke + antihyp + sbp + bio + diabetess + cutRt + offset(logw), data=ncc_s), 
warning=function(w) {tag <<-1})

if(tag==0) {GOF <- wald.test(vcov(NCC_s2),coef(NCC_s2), Terms=7:10)
GOF_ncc_s_wb <- cbind(as.vector(GOF$result$chi2)[1],"GOF_ncc_s_wb")
colnames(GOF_ncc_s_wb) <- c("GOF","name")} else {GOF_ncc_s_wb <- c(NA,"GOF_ncc_s_wb")}

  
#H-L
HL_ncc_s_wb <- cbind(HL(Uncc_s, Rt_ncc_s_wb_m,Uncc_s$CVD,w),"HL_ncc_s_wb")
colnames(HL_ncc_s_wb) <- c("HL","name")



#NRI 
NRI_ncc_s <- cbind(NRI(wc_s, wc_s$CVD, Rt_ncc_s_wb_m,Rt_ncc_s_wob_m, Uncc_s$CVD, c(0,0.05,0.20,100), w[Uncc_s$CVD==0]),"NRI_ncc_s")
colnames(NRI_ncc_s) <- c("NRI","name")




#IDI 
IDI_ncc_s <- cbind(IDI(Rt_ncc_s_wb_m,Rt_ncc_s_wob_m,Uncc_s$CVD,w), "IDI_ncc_s")
colnames(IDI_ncc_s) <- c("IDI","name")



#C-INDEX 
CIND_ncc_s_wb<-cbind(CIND(Rt_ncc_s_wb_m,Uncc_s$CVD,Uncc_s$stop,Uncc_s$TWINNR,w),"CIND_ncc_s_wb")
CIND_ncc_s_wob<-cbind(CIND(Rt_ncc_s_wob_m,Uncc_s$CVD,Uncc_s$stop,Uncc_s$TWINNR,w),"CIND_ncc_s_wob")
colnames(CIND_ncc_s_wb) <- c("cindex","name")
colnames(CIND_ncc_s_wob) <- c("cindex","name")


      
#################
#MERGING RESULTS#
#################

### a.ASSOCIATION ###

CO <- cbind(coWC, coPR,  coSP, coBR, coNCC, coBI_s4, coBII_s4, coCA_s4,coNCC_s,coNCC_s_br_age,coNCC_s_br_sex)
SE <- cbind(seWC, sePR,  seSP, seBR, seNCC, seBI_s4, seBII_s4, seCA_s4,seNCC_s)

#NUmber of non-repeated subjects
DIM <- round(cbind(L_random, L_random_ncc,  L_strata4, L_strata_ncc),0)


### b.INDIVIDUAL RISK ###

pred_Pr <- merge(Rt_wc,Rt_Pr, by="TWINNR")
pred_Pr_wc <- merge(Rt_wc,Rt_Pr_wc, by="TWINNR")
pred_ncc <- merge(Rt_wc,Rt_ncc, by="TWINNR")
pred_ncc_wc <- merge(Rt_wc,Rt_ncc_wc, by="TWINNR")
pred_BII <- merge(Rt_wc,Rt_BII, by="TWINNR")
pred_BII_wc <- merge(Rt_wc,Rt_BII_wc, by="TWINNR")
pred_ncc_s <- merge(Rt_wc,Rt_ncc_s, by="TWINNR")
pred_ncc_s_wc <- merge(Rt_wc,Rt_ncc_s_wc, by="TWINNR") 
IND <- list(pred_Pr,pred_Pr_wc,pred_ncc,pred_ncc_wc,pred_BII,pred_BII_wc,pred_ncc_s,pred_ncc_s_wc)

### c.PREDICTION MEASURES ###

GOFm <- rbind (GOF_wc_wb,GOF_Pr_wb,GOF_ncc_wb,GOF_BII_wb,GOF_ncc_s_wb)
HLm <- rbind(HL_wc_wb,HL_Pr_wb,HL_ncc_wb,HL_BII_wb,HL_ncc_s_wb)
NRIm <- rbind(NRI_wc,NRI_Pr,NRI_ncc,NRI_BII,NRI_ncc_s)
IDIm <- rbind(IDI_wc,IDI_Pr,IDI_ncc,IDI_BII,IDI_ncc_s)
CINDm <- rbind(CIND_wc_wb,CIND_wc_wob,CIND_Pr_wb,CIND_Pr_wob,CIND_ncc_wb,CIND_ncc_wob, CIND_BII_wb,CIND_BII_wob,CIND_ncc_s_wb,CIND_ncc_s_wob)


return(list(CO,SE, DIM, IND,GOFm,HLm,NRIm,IDIm,CINDm))


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
 

for (k in 1:nrep) {   F11<- ccfun(229,2,1095)
	                F11_beta <-data.frame(rbind(F11_beta, F11[[1]]))
	                F11_se <-data.frame(rbind(F11_se, F11[[2]]))
	                F11_size <-data.frame(rbind(F11_size, F11[[3]]))
	                 F11_Pr <-rbind(F11_Pr, F11[[4]][[1]])
	                F11_Pr_wc <-rbind(F11_Pr_wc, F11[[4]][[2]])
	                F11_ncc <-rbind(F11_ncc, F11[[4]][[3]])
	                F11_ncc_wc <-rbind(F11_ncc_wc, F11[[4]][[4]])
	                F11_BII <-rbind(F11_BII, F11[[4]][[5]])
	                F11_BII_wc <-rbind(F11_BII_wc, F11[[4]][[6]])
	                F11_ncc_s <-rbind(F11_ncc_s, F11[[4]][[7]])
	                F11_ncc_s_wc <-rbind(F11_ncc_s_wc, F11[[4]][[8]])
	                F11_GOF <- rbind(F11_GOF, F11[[5]])
	                F11_HL <-rbind(F11_HL, F11[[6]])
	                F11_NRI <-rbind(F11_NRI, F11[[7]])
	                F11_IDI <-rbind(F11_IDI, F11[[8]])
	                F11_CIND <-rbind(F11_CIND, F11[[9]])
                      print(k)
                      print(system("date"))}



save <- list(F11_beta,F11_se,F11_size,F11_Pr, F11_Pr_wc,F11_ncc,F11_ncc_wc, F11_BII,F11_BII_wc, F11_ncc_s, F11_ncc_s_wc, F11_GOF, F11_HL,F11_NRI,F11_IDI,F11_CIND)
save(save,file="Results/11_bio1.Rdata")


for (k in 1:nrep) {   F13<- ccfun(632,4,1095)
	                F13_beta <-data.frame(rbind(F13_beta, F13[[1]]))
	                F13_se <-data.frame(rbind(F13_se, F13[[2]]))
	                F13_size <-data.frame(rbind(F13_size, F13[[3]]))
	                 F13_Pr <-rbind(F13_Pr, F13[[4]][[1]])
	                F13_Pr_wc <-rbind(F13_Pr_wc, F13[[4]][[2]])
	                F13_ncc <-rbind(F13_ncc, F13[[4]][[3]])
	                F13_ncc_wc <-rbind(F13_ncc_wc, F13[[4]][[4]])
	                F13_BII <-rbind(F13_BII, F13[[4]][[5]])
	                F13_BII_wc <-rbind(F13_BII_wc, F13[[4]][[6]])
	                F13_ncc_s <-rbind(F13_ncc_s, F13[[4]][[7]])
	                F13_ncc_s_wc <-rbind(F13_ncc_s_wc, F13[[4]][[8]])
	                F13_GOF <- rbind(F13_GOF, F13[[5]])
	                F13_HL <-rbind(F13_HL, F13[[6]])
	                F13_NRI <-rbind(F13_NRI, F13[[7]])
	                F13_IDI <-rbind(F13_IDI, F13[[8]])
	                F13_CIND <-rbind(F13_CIND, F13[[9]])
                      print(k)
                      print(system("date"))}



save <- list(F13_beta,F13_se,F13_size,F13_Pr, F13_Pr_wc,F13_ncc,F13_ncc_wc, F13_BII,F13_BII_wc, F13_ncc_s, F13_ncc_s_wc, F13_GOF, F13_HL,F13_NRI,F13_IDI,F13_CIND)
save(save,file="Results/13_bio1.Rdata")

