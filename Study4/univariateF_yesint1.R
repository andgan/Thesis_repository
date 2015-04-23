#----------------------------------------------
# Filename: UnivariteF_yesint1.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in females. It obtaines association results for age-stratified analysis.
# Note: 
#-----------------------------------------------
# Data used: bdE9F.Rdata, AF1-AF5.Rdata
# Data created: Formulas: XFFbYI.Rdata,XFFCAYI.Rdata,XFFCVDYI.Rdata,XFFREYI.Rdata,XFFDGYI.Rdata,XFFECYI.Rdata,XFFODYI.Rdata,
#								Association results: XRESbYI.Rdata,XRESCAYI.Rdata,XRESCVDYI.Rdata,XRESREYI.Rdata,XRESDGYI.Rdata,XRESECYI.Rdata,XRESODYI.Rdata,XRESbHYI.Rdata
#								Number of deaths per category: XTABbYI.Rdata,XTABCAYI.Rdata,XTABCVDYI.Rdata,XTABREYI.Rdata,XTABDGYI.Rdata,XTABECYI.Rdata,XTABODYI.Rdata,XTABbHYI.Rdata
#-----------------------------------------------
# OP: R 3.1.2
#-----------------------------------------------*/



## Load packages
library(foreign)
library(survival)
library(rms)
library(mice)
library(glmnet)
library(lattice)
library(WGCNA)
library(Hmisc)
library(parallel)
library(doMC)
library(pec)
library(riskRegression)
library(cmprsk)

## Load in-house functions
source("/proj/b2011036/uk.biobank/Pgm/new_function.R")

load("/proj/b2011036/uk.biobank/imputation_results/bdE9F.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

N_MISSF <- bdE9F$nmis


#quantile(A1$age,seq(0,1,0.3333))
cutpoint <- c(53,62)

## CREATE AGE-SPLITTED DATASET TO TEST LACK OF HAZARD PROPORTIONALITY ##
AF_spl1 <- survSplit(AF1,cut=cutpoint,end="age_exit",event="out",start="age")
AF_spl2 <- survSplit(AF2,cut=cutpoint,end="age_exit",event="out",start="age")
AF_spl3 <- survSplit(AF3,cut=cutpoint,end="age_exit",event="out",start="age")
AF_spl4 <- survSplit(AF4,cut=cutpoint,end="age_exit",event="out",start="age")
AF_spl5 <- survSplit(AF5,cut=cutpoint,end="age_exit",event="out",start="age")

## CREAFE HEAVISIDE FUNCTIONS FOR TESTING LACK OF PROPORTIONALITY ##
AF_spl1$hv1=ifelse(AF_spl1$age<cutpoint[1],1,0)
AF_spl1$hv2=ifelse(AF_spl1$age<cutpoint[2] & AF_spl1$age>=cutpoint[1],1,0)
AF_spl1$hv3=ifelse(AF_spl1$age>=cutpoint[2],1,0)

AF_spl2$hv1=ifelse(AF_spl2$age<cutpoint[1],1,0)
AF_spl2$hv2=ifelse(AF_spl2$age<cutpoint[2] & AF_spl2$age>=cutpoint[1],1,0)
AF_spl2$hv3=ifelse(AF_spl2$age>=cutpoint[2],1,0)

AF_spl3$hv1=ifelse(AF_spl3$age<cutpoint[1],1,0)
AF_spl3$hv2=ifelse(AF_spl3$age<cutpoint[2] & AF_spl3$age>=cutpoint[1],1,0)
AF_spl3$hv3=ifelse(AF_spl3$age>=cutpoint[2],1,0)

AF_spl4$hv1=ifelse(AF_spl4$age<cutpoint[1],1,0)
AF_spl4$hv2=ifelse(AF_spl4$age<cutpoint[2] & AF_spl4$age>=cutpoint[1],1,0)
AF_spl4$hv3=ifelse(AF_spl4$age>=cutpoint[2],1,0)

AF_spl5$hv1=ifelse(AF_spl5$age<cutpoint[1],1,0)
AF_spl5$hv2=ifelse(AF_spl5$age<cutpoint[2] & AF_spl5$age>=cutpoint[1],1,0)
AF_spl5$hv3=ifelse(AF_spl5$age>=cutpoint[2],1,0)


## CAUSE-SPECIFIC MORTALITY ##
T1 <- AF1
T2 <- AF2
T3 <- AF3
T4 <- AF4
T5 <- AF5
nam <- c("CA","CVD","RE","DG","EC","OD")

for (l in nam)

{
	## CREATE AGE-SPLITTED DATASET TO TEST LACK OF HAZARD PROPORTIONALITY ##
	T1$event <- ifelse(AF1$status== (which(nam==l)),1,0)
	T2$event <- ifelse(AF2$status== (which(nam==l)),1,0)
	T3$event <- ifelse(AF3$status== (which(nam==l)),1,0)
	T4$event <- ifelse(AF4$status== (which(nam==l)),1,0)
	T5$event <- ifelse(AF5$status== (which(nam==l)),1,0)

	T_spl1 <- survSplit(T1,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl2 <- survSplit(T2,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl3 <- survSplit(T3,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl4 <- survSplit(T4,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl5 <- survSplit(T5,cut=cutpoint,end="age_exit",event="event",start="age")

	## CREAFE HEAVISIDE FUNCTIONS FOR TESTING LACK OF PROPORTIONALITY ##
	T_spl1$hv1=ifelse(T_spl1$age<cutpoint[1],1,0)
	T_spl1$hv2=ifelse(T_spl1$age<cutpoint[2] & T_spl1$age>=cutpoint[1],1,0)
	T_spl1$hv3=ifelse(T_spl1$age>=cutpoint[2],1,0)

	T_spl2$hv1=ifelse(T_spl2$age<cutpoint[1],1,0)
	T_spl2$hv2=ifelse(T_spl2$age<cutpoint[2] & T_spl2$age>=cutpoint[1],1,0)
	T_spl2$hv3=ifelse(T_spl2$age>=cutpoint[2],1,0)

	T_spl3$hv1=ifelse(T_spl3$age<cutpoint[1],1,0)
	T_spl3$hv2=ifelse(T_spl3$age<cutpoint[2] & T_spl3$age>=cutpoint[1],1,0)
	T_spl3$hv3=ifelse(T_spl3$age>=cutpoint[2],1,0)

	T_spl4$hv1=ifelse(T_spl4$age<cutpoint[1],1,0)
	T_spl4$hv2=ifelse(T_spl4$age<cutpoint[2] & T_spl4$age>=cutpoint[1],1,0)
	T_spl4$hv3=ifelse(T_spl4$age>=cutpoint[2],1,0)

	T_spl5$hv1=ifelse(T_spl5$age<cutpoint[1],1,0)
	T_spl5$hv2=ifelse(T_spl5$age<cutpoint[2] & T_spl5$age>=cutpoint[1],1,0)
	T_spl5$hv3=ifelse(T_spl5$age>=cutpoint[2],1,0)

	# Assign T_spl'x' to the right dataset
	assign(paste0("AF",l,"_spl1"), T_spl1)
	assign(paste0("AF",l,"_spl2"), T_spl2)
	assign(paste0("AF",l,"_spl3"), T_spl3)
	assign(paste0("AF",l,"_spl4"), T_spl4)
	assign(paste0("AF",l,"_spl5"), T_spl5)
}

rm(T1,T2,T3,T4,T5)

## CREATE AN EXTERNAL VALIDATION SET ##
set.seed(123)
idsplitF2 <- which(AF1$center%in%c("Glasgow","Edinburgh"))

AFT1 <- AF1[-unique(c(idsplitF2)),]
AFT2 <- AF2[-unique(c(idsplitF2)),]
AFT3 <- AF3[-unique(c(idsplitF2)),]
AFT4 <- AF4[-unique(c(idsplitF2)),]
AFT5 <- AF5[-unique(c(idsplitF2)),]


# Name of variables to include in the analysis
NAMESF <- colnames(AF1)[!colnames(AF1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","charlsonSR")]


# Set variables
XFFb <- NULL
XFFCA <- NULL
XFFCVD <- NULL
XFFRE <- NULL
XFFDG <- NULL
XFFEC <- NULL
XFFOD <- NULL

R2F <- NULL
PROPpF <- NULL

XRESb <- list()
XRESCA <- list()
XRESCVD <- list()
XRESRE <- list()
XRESDG <- list()
XRESEC <- list()
XRESOD <- list()
XRESbH <- list()

XTABb <- list()
XTABCA <- list()
XTABCVD <- list()
XTABRE <- list()
XTABDG <- list()
XTABEC <- list()
XTABOD <- list()
XTABbH <- list()

##########################
#### 3b. RUN ANALYSIS ####
##########################

for (k in NAMESF)
{	

	# Set variables to -999, which will be assigned if not any other value is assigned

	FFb <- -999
	FFCA <- -999
	FFCVD <- -999
	FFRE <- -999
	FFDG <- -999
	FFEC <- -999
	FFOD <- -999

	RESb <- -999
	RESCA <- -999
	RESCVD <- -999
	RESRE <- -999
	RESDG <- -999
	RESEC <- -999
	RESOD <- -999
	RESbH <- -999

	
	# Set reference level to most frequent level
	tt <- table(AF1[,k])
	AF_spl1[,k]  <- relevel(AF_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AF_spl2[,k]  <- relevel(AF_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AF_spl3[,k]  <- relevel(AF_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AF_spl4[,k]  <- relevel(AF_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AF_spl5[,k]  <- relevel(AF_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	# Add a new variable which is called hv'i'_'k'
	for (i in 1:3)
	{
		AF_spl1[paste0("hv",i,"_",k)] <- as.factor(AF_spl1[,paste0("hv",i)]*(as.numeric(AF_spl1[,k])-1))
		AF_spl2[paste0("hv",i,"_",k)] <- as.factor(AF_spl2[,paste0("hv",i)]*(as.numeric(AF_spl2[,k])-1))
		AF_spl3[paste0("hv",i,"_",k)] <- as.factor(AF_spl3[,paste0("hv",i)]*(as.numeric(AF_spl3[,k])-1))
		AF_spl4[paste0("hv",i,"_",k)] <- as.factor(AF_spl4[,paste0("hv",i)]*(as.numeric(AF_spl4[,k])-1))
		AF_spl5[paste0("hv",i,"_",k)] <- as.factor(AF_spl5[,paste0("hv",i)]*(as.numeric(AF_spl5[,k])-1))			
	}

  # Next it sets the reference level to most frequent level for all the cause-specific mortality
  # And add a new variable which is called hv'i'_'k'
	AFCA_spl1[,k]  <- relevel(AFCA_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFCA_spl2[,k]  <- relevel(AFCA_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFCA_spl3[,k]  <- relevel(AFCA_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFCA_spl4[,k]  <- relevel(AFCA_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFCA_spl5[,k]  <- relevel(AFCA_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFCA_spl1[paste0("hv",i,"_",k)] <- as.factor(AFCA_spl1[,paste0("hv",i)]*(as.numeric(AFCA_spl1[,k])-1))
		AFCA_spl2[paste0("hv",i,"_",k)] <- as.factor(AFCA_spl2[,paste0("hv",i)]*(as.numeric(AFCA_spl2[,k])-1))
		AFCA_spl3[paste0("hv",i,"_",k)] <- as.factor(AFCA_spl3[,paste0("hv",i)]*(as.numeric(AFCA_spl3[,k])-1))
		AFCA_spl4[paste0("hv",i,"_",k)] <- as.factor(AFCA_spl4[,paste0("hv",i)]*(as.numeric(AFCA_spl4[,k])-1))
		AFCA_spl5[paste0("hv",i,"_",k)] <- as.factor(AFCA_spl5[,paste0("hv",i)]*(as.numeric(AFCA_spl5[,k])-1))			
	}

	AFCVD_spl1[,k]  <- relevel(AFCVD_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFCVD_spl2[,k]  <- relevel(AFCVD_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFCVD_spl3[,k]  <- relevel(AFCVD_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFCVD_spl4[,k]  <- relevel(AFCVD_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFCVD_spl5[,k]  <- relevel(AFCVD_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFCVD_spl1[paste0("hv",i,"_",k)] <- as.factor(AFCVD_spl1[,paste0("hv",i)]*(as.numeric(AFCVD_spl1[,k])-1))
		AFCVD_spl2[paste0("hv",i,"_",k)] <- as.factor(AFCVD_spl2[,paste0("hv",i)]*(as.numeric(AFCVD_spl2[,k])-1))
		AFCVD_spl3[paste0("hv",i,"_",k)] <- as.factor(AFCVD_spl3[,paste0("hv",i)]*(as.numeric(AFCVD_spl3[,k])-1))
		AFCVD_spl4[paste0("hv",i,"_",k)] <- as.factor(AFCVD_spl4[,paste0("hv",i)]*(as.numeric(AFCVD_spl4[,k])-1))
		AFCVD_spl5[paste0("hv",i,"_",k)] <- as.factor(AFCVD_spl5[,paste0("hv",i)]*(as.numeric(AFCVD_spl5[,k])-1))			
	}

	AFRE_spl1[,k]  <- relevel(AFRE_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFRE_spl2[,k]  <- relevel(AFRE_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFRE_spl3[,k]  <- relevel(AFRE_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFRE_spl4[,k]  <- relevel(AFRE_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFRE_spl5[,k]  <- relevel(AFRE_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFRE_spl1[paste0("hv",i,"_",k)] <- as.factor(AFRE_spl1[,paste0("hv",i)]*(as.numeric(AFRE_spl1[,k])-1))
		AFRE_spl2[paste0("hv",i,"_",k)] <- as.factor(AFRE_spl2[,paste0("hv",i)]*(as.numeric(AFRE_spl2[,k])-1))
		AFRE_spl3[paste0("hv",i,"_",k)] <- as.factor(AFRE_spl3[,paste0("hv",i)]*(as.numeric(AFRE_spl3[,k])-1))
		AFRE_spl4[paste0("hv",i,"_",k)] <- as.factor(AFRE_spl4[,paste0("hv",i)]*(as.numeric(AFRE_spl4[,k])-1))
		AFRE_spl5[paste0("hv",i,"_",k)] <- as.factor(AFRE_spl5[,paste0("hv",i)]*(as.numeric(AFRE_spl5[,k])-1))			
	}

	AFDG_spl1[,k]  <- relevel(AFDG_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFDG_spl2[,k]  <- relevel(AFDG_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFDG_spl3[,k]  <- relevel(AFDG_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFDG_spl4[,k]  <- relevel(AFDG_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFDG_spl5[,k]  <- relevel(AFDG_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFDG_spl1[paste0("hv",i,"_",k)] <- as.factor(AFDG_spl1[,paste0("hv",i)]*(as.numeric(AFDG_spl1[,k])-1))
		AFDG_spl2[paste0("hv",i,"_",k)] <- as.factor(AFDG_spl2[,paste0("hv",i)]*(as.numeric(AFDG_spl2[,k])-1))
		AFDG_spl3[paste0("hv",i,"_",k)] <- as.factor(AFDG_spl3[,paste0("hv",i)]*(as.numeric(AFDG_spl3[,k])-1))
		AFDG_spl4[paste0("hv",i,"_",k)] <- as.factor(AFDG_spl4[,paste0("hv",i)]*(as.numeric(AFDG_spl4[,k])-1))
		AFDG_spl5[paste0("hv",i,"_",k)] <- as.factor(AFDG_spl5[,paste0("hv",i)]*(as.numeric(AFDG_spl5[,k])-1))			
	}

	AFEC_spl1[,k]  <- relevel(AFEC_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFEC_spl2[,k]  <- relevel(AFEC_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFEC_spl3[,k]  <- relevel(AFEC_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFEC_spl4[,k]  <- relevel(AFEC_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFEC_spl5[,k]  <- relevel(AFEC_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFEC_spl1[paste0("hv",i,"_",k)] <- as.factor(AFEC_spl1[,paste0("hv",i)]*(as.numeric(AFEC_spl1[,k])-1))
		AFEC_spl2[paste0("hv",i,"_",k)] <- as.factor(AFEC_spl2[,paste0("hv",i)]*(as.numeric(AFEC_spl2[,k])-1))
		AFEC_spl3[paste0("hv",i,"_",k)] <- as.factor(AFEC_spl3[,paste0("hv",i)]*(as.numeric(AFEC_spl3[,k])-1))
		AFEC_spl4[paste0("hv",i,"_",k)] <- as.factor(AFEC_spl4[,paste0("hv",i)]*(as.numeric(AFEC_spl4[,k])-1))
		AFEC_spl5[paste0("hv",i,"_",k)] <- as.factor(AFEC_spl5[,paste0("hv",i)]*(as.numeric(AFEC_spl5[,k])-1))			
	}

	AFOD_spl1[,k]  <- relevel(AFOD_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AFOD_spl2[,k]  <- relevel(AFOD_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AFOD_spl3[,k]  <- relevel(AFOD_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AFOD_spl4[,k]  <- relevel(AFOD_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AFOD_spl5[,k]  <- relevel(AFOD_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AFOD_spl1[paste0("hv",i,"_",k)] <- as.factor(AFOD_spl1[,paste0("hv",i)]*(as.numeric(AFOD_spl1[,k])-1))
		AFOD_spl2[paste0("hv",i,"_",k)] <- as.factor(AFOD_spl2[,paste0("hv",i)]*(as.numeric(AFOD_spl2[,k])-1))
		AFOD_spl3[paste0("hv",i,"_",k)] <- as.factor(AFOD_spl3[,paste0("hv",i)]*(as.numeric(AFOD_spl3[,k])-1))
		AFOD_spl4[paste0("hv",i,"_",k)] <- as.factor(AFOD_spl4[,paste0("hv",i)]*(as.numeric(AFOD_spl4[,k])-1))
		AFOD_spl5[paste0("hv",i,"_",k)] <- as.factor(AFOD_spl5[,paste0("hv",i)]*(as.numeric(AFOD_spl5[,k])-1))			
	}

	# N.of events in each class #
	tab1 <- table(eval(parse(text=paste0("AF_spl1$","hv1_",k,"[AF_spl1$out==1]"))))
	tab2 <- table(eval(parse(text=paste0("AF_spl1$","hv2_",k,"[AF_spl1$out==1]"))))
	tab3 <- table(eval(parse(text=paste0("AF_spl1$","hv3_",k,"[AF_spl1$out==1]"))))
	tabb <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFCA_spl1$","hv1_",k,"[AFCA_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFCA_spl1$","hv2_",k,"[AFCA_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFCA_spl1$","hv3_",k,"[AFCA_spl1$event==1]"))))
	tabCA <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFCVD_spl1$","hv1_",k,"[AFCVD_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFCVD_spl1$","hv2_",k,"[AFCVD_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFCVD_spl1$","hv3_",k,"[AFCVD_spl1$event==1]"))))
	tabCVD <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFRE_spl1$","hv1_",k,"[AFRE_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFRE_spl1$","hv2_",k,"[AFRE_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFRE_spl1$","hv3_",k,"[AFRE_spl1$event==1]"))))
	tabRE <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFDG_spl1$","hv1_",k,"[AFDG_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFDG_spl1$","hv2_",k,"[AFDG_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFDG_spl1$","hv3_",k,"[AFDG_spl1$event==1]"))))
	tabDG <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFEC_spl1$","hv1_",k,"[AFEC_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFEC_spl1$","hv2_",k,"[AFEC_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFEC_spl1$","hv3_",k,"[AFEC_spl1$event==1]"))))
	tabEC <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AFOD_spl1$","hv1_",k,"[AFOD_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AFOD_spl1$","hv2_",k,"[AFOD_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AFOD_spl1$","hv3_",k,"[AFOD_spl1$event==1]"))))
	tabOD <- rbind(tab1,tab2,tab3)
	

	# Fodel for overall mortality - time-dependent #
	modTIb <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AF_spl"), 
		error=function(e) {NULL})

		AFH_spl1 <- AF_spl1[AF_spl1$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
		AFH_spl2 <- AF_spl2[AF_spl2$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
		AFH_spl3 <- AF_spl3[AF_spl3$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
		AFH_spl4 <- AF_spl4[AF_spl4$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
		AFH_spl5 <- AF_spl5[AF_spl5$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]

		tab1 <- table(eval(parse(text=paste0("AFH_spl1$","hv1_",k,"[AFH_spl1$out==1]"))))
		tab2 <- table(eval(parse(text=paste0("AFH_spl1$","hv2_",k,"[AFH_spl1$out==1]"))))
		tab3 <- table(eval(parse(text=paste0("AFH_spl1$","hv3_",k,"[AFH_spl1$out==1]"))))
		tabbH <- rbind(tab1,tab2,tab3)


		modTIbH <- tryCatch(
			coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFH_spl"), 
			error=function(e) {NULL})


	# Fodels for cause-specific mortality - time-dependent #
	modTICA <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFCA_spl"), 
		error=function(e) {NULL})

	modTICVD <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFCVD_spl"), 
		error=function(e) {NULL})

	modTIRE <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFRE_spl"), 
		error=function(e) {NULL})

	modTIDG <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFDG_spl"), 
		error=function(e) {NULL})

	modTIEC <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFEC_spl"), 
		error=function(e) {NULL})

	modTIOD <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSF, "AFOD_spl"), 
		error=function(e) {NULL})	

			
	if(!is.null(modTIb))
	{
		FFb <- paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <-  tryCatch(pool(modTIb), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESb <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		} 			
	}	

	if(!is.null(modTIbH))
	{
		FFb <- paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <-  tryCatch(pool(modTIbH), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESbH <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		} 				
	}
	
	if(!is.null(modTICA))
	{
		FFCA <- paste0("Surv(age,age_exit,status==1)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTICA), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESCA <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}						
	}

	if(!is.null(modTICVD))
	{
		FFCVD <- paste0("Surv(age,age_exit,status==2)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTICVD), error=function(e) {NULL})	
		if (!is.null(pp))
		{
			RESCVD <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]), c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}				
	}

	if(!is.null(modTIRE))
	{
		FFRE <- paste0("Surv(age,age_exit,status==3)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIRE), error=function(e) {NULL})	
		if (!is.null(pp))
		{
			RESRE <- cbind(rep(levels(AF_spl1[,k]),3), c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}			
	}

	if(!is.null(modTIDG))
	{
		FFDG <- paste0("Surv(age,age_exit,status==4)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIDG), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESDG <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}				
	}

	if(!is.null(modTIEC))
	{
		FFEC <- paste0("Surv(age,age_exit,status==5)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIEC), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESEC <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}			
	}

	if(!is.null(modTIOD))
	{
		FFOD <- paste0("Surv(age,age_exit,status==6)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIOD), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESOD <- cbind(rep(levels(AF_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]), c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))	
		}						
	}	

	##  Collect the results ##
	
	# Formulas
	XFFb <- rbind(XFFb,FFb)
	XFFCA <- rbind(XFFCA,FFCA)
	XFFCVD <- rbind(XFFCVD,FFCVD)
	XFFRE <- rbind(XFFRE,FFRE)
	XFFDG <- rbind(XFFDG,FFDG)
	XFFEC <- rbind(XFFEC,FFEC)
	XFFOD <- rbind(XFFOD,FFOD)

	# Coefficients + standard errors
	XRESb[[length(XRESb)+1]] <- RESb
	XRESCA[[length(XRESCA)+1]] <- RESCA
	XRESCVD[[length(XRESCVD)+1]] <- RESCVD
	XRESRE[[length(XRESRE)+1]] <- RESRE
	XRESDG[[length(XRESDG)+1]] <- RESDG
	XRESEC[[length(XRESEC)+1]] <- RESEC
	XRESOD[[length(XRESOD)+1]] <- RESOD
	
	XRESbH[[length(XRESbH)+1]] <- RESbH
	
	XTABb[[length(XTABb)+1]] <- tabb
	XTABCA[[length(XTABCA)+1]] <- tabCA
	XTABCVD[[length(XTABCVD)+1]] <- tabCVD
	XTABRE[[length(XTABRE)+1]] <- tabRE
	XTABDG[[length(XTABDG)+1]] <- tabDG
	XTABEC[[length(XTABEC)+1]] <- tabEC
	XTABOD[[length(XTABOD)+1]] <- tabOD
	
	XTABbH[[length(XTABbH)+1]] <- tabbH
	
	print(k)
}

## SAVE OBJECTS ##
for (i in c("XFFb","XFFCA","XFFCVD","XFFRE","XFFDG","XFFEC","XFFOD",
"XRESb","XRESCA","XRESCVD","XRESRE","XRESDG","XRESEC","XRESOD","XRESbH",
"XTABb","XTABCA","XTABCVD","XTABRE","XTABDG","XTABEC","XTABOD","XTABbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))}