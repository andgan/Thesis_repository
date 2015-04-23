#----------------------------------------------
# Filename: UnivariteM_yesint1.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in males. It obtaines association results for age-stratified analysis.
# Note: 
#-----------------------------------------------
# Data used: bdE9M.Rdata, AM1-AM5.Rdata
# Data created: Formulas: YFFbYI.Rdata,YFFCAYI.Rdata,YFFCVDYI.Rdata,YFFREYI.Rdata,YFFDGYI.Rdata,YFFECYI.Rdata,YFFODYI.Rdata,
#								Association results: YRESbYI.Rdata,YRESCAYI.Rdata,YRESCVDYI.Rdata,YRESREYI.Rdata,YRESDGYI.Rdata,YRESECYI.Rdata,YRESODYI.Rdata,YRESbHYI.Rdata
#								Number of deaths per category: YTABbYI.Rdata,YTABCAYI.Rdata,YTABCVDYI.Rdata,YTABREYI.Rdata,YTABDGYI.Rdata,YTABECYI.Rdata,YTABODYI.Rdata,YTABbHYI.Rdata
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


load("/proj/b2011036/uk.biobank/imputation_results/bdE9M.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")

N_MISSM <- bdE9M$nmis



#quantile(AM1$age,seq(0,1,0.3333))
cutpoint <- c(53,62)

## CREATE AGE-SPLITTED DATASET TO TEST LACK OF HAZARD PROPORTIONALITY ##
AM_spl1 <- survSplit(AM1,cut=cutpoint,end="age_exit",event="out",start="age")
AM_spl2 <- survSplit(AM2,cut=cutpoint,end="age_exit",event="out",start="age")
AM_spl3 <- survSplit(AM3,cut=cutpoint,end="age_exit",event="out",start="age")
AM_spl4 <- survSplit(AM4,cut=cutpoint,end="age_exit",event="out",start="age")
AM_spl5 <- survSplit(AM5,cut=cutpoint,end="age_exit",event="out",start="age")

## CREAME HEAVISIDE FUNCTIONS FOR TESTING LACK OF PROPORTIONALITY ##
AM_spl1$hv1=ifelse(AM_spl1$age<cutpoint[1],1,0)
AM_spl1$hv2=ifelse(AM_spl1$age<cutpoint[2] & AM_spl1$age>=cutpoint[1],1,0)
AM_spl1$hv3=ifelse(AM_spl1$age>=cutpoint[2],1,0)

AM_spl2$hv1=ifelse(AM_spl2$age<cutpoint[1],1,0)
AM_spl2$hv2=ifelse(AM_spl2$age<cutpoint[2] & AM_spl2$age>=cutpoint[1],1,0)
AM_spl2$hv3=ifelse(AM_spl2$age>=cutpoint[2],1,0)

AM_spl3$hv1=ifelse(AM_spl3$age<cutpoint[1],1,0)
AM_spl3$hv2=ifelse(AM_spl3$age<cutpoint[2] & AM_spl3$age>=cutpoint[1],1,0)
AM_spl3$hv3=ifelse(AM_spl3$age>=cutpoint[2],1,0)

AM_spl4$hv1=ifelse(AM_spl4$age<cutpoint[1],1,0)
AM_spl4$hv2=ifelse(AM_spl4$age<cutpoint[2] & AM_spl4$age>=cutpoint[1],1,0)
AM_spl4$hv3=ifelse(AM_spl4$age>=cutpoint[2],1,0)

AM_spl5$hv1=ifelse(AM_spl5$age<cutpoint[1],1,0)
AM_spl5$hv2=ifelse(AM_spl5$age<cutpoint[2] & AM_spl5$age>=cutpoint[1],1,0)
AM_spl5$hv3=ifelse(AM_spl5$age>=cutpoint[2],1,0)


## CAUSE-SPECIFIC MORTALITY ##
T1 <- AM1
T2 <- AM2
T3 <- AM3
T4 <- AM4
T5 <- AM5
nam <- c("CA","CVD","RE","DG","EC","OD")

for (l in nam)

{
	## CREATE AGE-SPLITTED DATASET TO TEST LACK OF HAZARD PROPORTIONALITY ##
	T1$event <- ifelse(AM1$status== (which(nam==l)),1,0)
	T2$event <- ifelse(AM2$status== (which(nam==l)),1,0)
	T3$event <- ifelse(AM3$status== (which(nam==l)),1,0)
	T4$event <- ifelse(AM4$status== (which(nam==l)),1,0)
	T5$event <- ifelse(AM5$status== (which(nam==l)),1,0)

	T_spl1 <- survSplit(T1,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl2 <- survSplit(T2,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl3 <- survSplit(T3,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl4 <- survSplit(T4,cut=cutpoint,end="age_exit",event="event",start="age")
	T_spl5 <- survSplit(T5,cut=cutpoint,end="age_exit",event="event",start="age")

	## CREAME HEAVISIDE FUNCTIONS FOR TESTING LACK OF PROPORTIONALITY ##
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
	assign(paste0("AM",l,"_spl1"), T_spl1)
	assign(paste0("AM",l,"_spl2"), T_spl2)
	assign(paste0("AM",l,"_spl3"), T_spl3)
	assign(paste0("AM",l,"_spl4"), T_spl4)
	assign(paste0("AM",l,"_spl5"), T_spl5)
}


## CREATE AN EXTERNAL VALIDATION SET ##
set.seed(123)
idsplitM2 <- which(AM1$center%in%c("Glasgow","Edinburgh"))

AMT1 <- AM1[-unique(c(idsplitM2)),]
AMT2 <- AM2[-unique(c(idsplitM2)),]
AMT3 <- AM3[-unique(c(idsplitM2)),]
AMT4 <- AM4[-unique(c(idsplitM2)),]
AMT5 <- AM5[-unique(c(idsplitM2)),]



# Name of variables to include in the analysis
NAMESM <- colnames(AM1)[!colnames(AM1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","charlsonSR")]


# Set variables
YFFb <- NULL
YFFCA <- NULL
YFFCVD <- NULL
YFFRE <- NULL
YFFDG <- NULL
YFFEC <- NULL
YFFOD <- NULL

R2M <- NULL
PROPpM <- NULL

YRESb <- list()
YRESCA <- list()
YRESCVD <- list()
YRESRE <- list()
YRESDG <- list()
YRESEC <- list()
YRESOD <- list()
YRESbH <- list()

YTABb <- list()
YTABCA <- list()
YTABCVD <- list()
YTABRE <- list()
YTABDG <- list()
YTABEC <- list()
YTABOD <- list()
YTABbH <- list()

##########################
#### 2b. RUN ANALYSIS ####
##########################

for (k in NAMESM)
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
	
	ptm <- proc.time()
	
		
	# Set reference level to most frequent level
	tt <- table(AM1[,k])
	AM_spl1[,k]  <- relevel(AM_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AM_spl2[,k]  <- relevel(AM_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AM_spl3[,k]  <- relevel(AM_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AM_spl4[,k]  <- relevel(AM_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AM_spl5[,k]  <- relevel(AM_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	# Add a new variable which is called hv'i'_'k'
	for (i in 1:3)
	{
		AM_spl1[paste0("hv",i,"_",k)] <- as.factor(AM_spl1[,paste0("hv",i)]*(as.numeric(AM_spl1[,k])-1))
		AM_spl2[paste0("hv",i,"_",k)] <- as.factor(AM_spl2[,paste0("hv",i)]*(as.numeric(AM_spl2[,k])-1))
		AM_spl3[paste0("hv",i,"_",k)] <- as.factor(AM_spl3[,paste0("hv",i)]*(as.numeric(AM_spl3[,k])-1))
		AM_spl4[paste0("hv",i,"_",k)] <- as.factor(AM_spl4[,paste0("hv",i)]*(as.numeric(AM_spl4[,k])-1))
		AM_spl5[paste0("hv",i,"_",k)] <- as.factor(AM_spl5[,paste0("hv",i)]*(as.numeric(AM_spl5[,k])-1))			
	}

  # Next it sets the reference level to most frequent level for all the cause-specific mortality
  # And add a new variable which is called hv'i'_'k'
	AMCA_spl1[,k]  <- relevel(AMCA_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMCA_spl2[,k]  <- relevel(AMCA_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMCA_spl3[,k]  <- relevel(AMCA_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMCA_spl4[,k]  <- relevel(AMCA_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMCA_spl5[,k]  <- relevel(AMCA_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMCA_spl1[paste0("hv",i,"_",k)] <- as.factor(AMCA_spl1[,paste0("hv",i)]*(as.numeric(AMCA_spl1[,k])-1))
		AMCA_spl2[paste0("hv",i,"_",k)] <- as.factor(AMCA_spl2[,paste0("hv",i)]*(as.numeric(AMCA_spl2[,k])-1))
		AMCA_spl3[paste0("hv",i,"_",k)] <- as.factor(AMCA_spl3[,paste0("hv",i)]*(as.numeric(AMCA_spl3[,k])-1))
		AMCA_spl4[paste0("hv",i,"_",k)] <- as.factor(AMCA_spl4[,paste0("hv",i)]*(as.numeric(AMCA_spl4[,k])-1))
		AMCA_spl5[paste0("hv",i,"_",k)] <- as.factor(AMCA_spl5[,paste0("hv",i)]*(as.numeric(AMCA_spl5[,k])-1))			
	}

	AMCVD_spl1[,k]  <- relevel(AMCVD_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMCVD_spl2[,k]  <- relevel(AMCVD_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMCVD_spl3[,k]  <- relevel(AMCVD_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMCVD_spl4[,k]  <- relevel(AMCVD_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMCVD_spl5[,k]  <- relevel(AMCVD_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMCVD_spl1[paste0("hv",i,"_",k)] <- as.factor(AMCVD_spl1[,paste0("hv",i)]*(as.numeric(AMCVD_spl1[,k])-1))
		AMCVD_spl2[paste0("hv",i,"_",k)] <- as.factor(AMCVD_spl2[,paste0("hv",i)]*(as.numeric(AMCVD_spl2[,k])-1))
		AMCVD_spl3[paste0("hv",i,"_",k)] <- as.factor(AMCVD_spl3[,paste0("hv",i)]*(as.numeric(AMCVD_spl3[,k])-1))
		AMCVD_spl4[paste0("hv",i,"_",k)] <- as.factor(AMCVD_spl4[,paste0("hv",i)]*(as.numeric(AMCVD_spl4[,k])-1))
		AMCVD_spl5[paste0("hv",i,"_",k)] <- as.factor(AMCVD_spl5[,paste0("hv",i)]*(as.numeric(AMCVD_spl5[,k])-1))			
	}

	AMRE_spl1[,k]  <- relevel(AMRE_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMRE_spl2[,k]  <- relevel(AMRE_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMRE_spl3[,k]  <- relevel(AMRE_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMRE_spl4[,k]  <- relevel(AMRE_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMRE_spl5[,k]  <- relevel(AMRE_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMRE_spl1[paste0("hv",i,"_",k)] <- as.factor(AMRE_spl1[,paste0("hv",i)]*(as.numeric(AMRE_spl1[,k])-1))
		AMRE_spl2[paste0("hv",i,"_",k)] <- as.factor(AMRE_spl2[,paste0("hv",i)]*(as.numeric(AMRE_spl2[,k])-1))
		AMRE_spl3[paste0("hv",i,"_",k)] <- as.factor(AMRE_spl3[,paste0("hv",i)]*(as.numeric(AMRE_spl3[,k])-1))
		AMRE_spl4[paste0("hv",i,"_",k)] <- as.factor(AMRE_spl4[,paste0("hv",i)]*(as.numeric(AMRE_spl4[,k])-1))
		AMRE_spl5[paste0("hv",i,"_",k)] <- as.factor(AMRE_spl5[,paste0("hv",i)]*(as.numeric(AMRE_spl5[,k])-1))			
	}

	AMDG_spl1[,k]  <- relevel(AMDG_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMDG_spl2[,k]  <- relevel(AMDG_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMDG_spl3[,k]  <- relevel(AMDG_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMDG_spl4[,k]  <- relevel(AMDG_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMDG_spl5[,k]  <- relevel(AMDG_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMDG_spl1[paste0("hv",i,"_",k)] <- as.factor(AMDG_spl1[,paste0("hv",i)]*(as.numeric(AMDG_spl1[,k])-1))
		AMDG_spl2[paste0("hv",i,"_",k)] <- as.factor(AMDG_spl2[,paste0("hv",i)]*(as.numeric(AMDG_spl2[,k])-1))
		AMDG_spl3[paste0("hv",i,"_",k)] <- as.factor(AMDG_spl3[,paste0("hv",i)]*(as.numeric(AMDG_spl3[,k])-1))
		AMDG_spl4[paste0("hv",i,"_",k)] <- as.factor(AMDG_spl4[,paste0("hv",i)]*(as.numeric(AMDG_spl4[,k])-1))
		AMDG_spl5[paste0("hv",i,"_",k)] <- as.factor(AMDG_spl5[,paste0("hv",i)]*(as.numeric(AMDG_spl5[,k])-1))			
	}

	AMEC_spl1[,k]  <- relevel(AMEC_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMEC_spl2[,k]  <- relevel(AMEC_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMEC_spl3[,k]  <- relevel(AMEC_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMEC_spl4[,k]  <- relevel(AMEC_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMEC_spl5[,k]  <- relevel(AMEC_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMEC_spl1[paste0("hv",i,"_",k)] <- as.factor(AMEC_spl1[,paste0("hv",i)]*(as.numeric(AMEC_spl1[,k])-1))
		AMEC_spl2[paste0("hv",i,"_",k)] <- as.factor(AMEC_spl2[,paste0("hv",i)]*(as.numeric(AMEC_spl2[,k])-1))
		AMEC_spl3[paste0("hv",i,"_",k)] <- as.factor(AMEC_spl3[,paste0("hv",i)]*(as.numeric(AMEC_spl3[,k])-1))
		AMEC_spl4[paste0("hv",i,"_",k)] <- as.factor(AMEC_spl4[,paste0("hv",i)]*(as.numeric(AMEC_spl4[,k])-1))
		AMEC_spl5[paste0("hv",i,"_",k)] <- as.factor(AMEC_spl5[,paste0("hv",i)]*(as.numeric(AMEC_spl5[,k])-1))			
	}

	AMOD_spl1[,k]  <- relevel(AMOD_spl1[,k], ref=names(tt)[tt==max(tt)][1])
	AMOD_spl2[,k]  <- relevel(AMOD_spl2[,k], ref=names(tt)[tt==max(tt)][1])
	AMOD_spl3[,k]  <- relevel(AMOD_spl3[,k], ref=names(tt)[tt==max(tt)][1])
	AMOD_spl4[,k]  <- relevel(AMOD_spl4[,k], ref=names(tt)[tt==max(tt)][1])
	AMOD_spl5[,k]  <- relevel(AMOD_spl5[,k], ref=names(tt)[tt==max(tt)][1])
	
	for (i in 1:3)
	{
		AMOD_spl1[paste0("hv",i,"_",k)] <- as.factor(AMOD_spl1[,paste0("hv",i)]*(as.numeric(AMOD_spl1[,k])-1))
		AMOD_spl2[paste0("hv",i,"_",k)] <- as.factor(AMOD_spl2[,paste0("hv",i)]*(as.numeric(AMOD_spl2[,k])-1))
		AMOD_spl3[paste0("hv",i,"_",k)] <- as.factor(AMOD_spl3[,paste0("hv",i)]*(as.numeric(AMOD_spl3[,k])-1))
		AMOD_spl4[paste0("hv",i,"_",k)] <- as.factor(AMOD_spl4[,paste0("hv",i)]*(as.numeric(AMOD_spl4[,k])-1))
		AMOD_spl5[paste0("hv",i,"_",k)] <- as.factor(AMOD_spl5[,paste0("hv",i)]*(as.numeric(AMOD_spl5[,k])-1))			
	}

	# N.of events in each class #		
	tab1 <- table(eval(parse(text=paste0("AM_spl1$","hv1_",k,"[AM_spl1$out==1]"))))
	tab2 <- table(eval(parse(text=paste0("AM_spl1$","hv2_",k,"[AM_spl1$out==1]"))))
	tab3 <- table(eval(parse(text=paste0("AM_spl1$","hv3_",k,"[AM_spl1$out==1]"))))
	tabb <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMCA_spl1$","hv1_",k,"[AMCA_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMCA_spl1$","hv2_",k,"[AMCA_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMCA_spl1$","hv3_",k,"[AMCA_spl1$event==1]"))))
	tabCA <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMCVD_spl1$","hv1_",k,"[AMCVD_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMCVD_spl1$","hv2_",k,"[AMCVD_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMCVD_spl1$","hv3_",k,"[AMCVD_spl1$event==1]"))))
	tabCVD <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMRE_spl1$","hv1_",k,"[AMRE_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMRE_spl1$","hv2_",k,"[AMRE_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMRE_spl1$","hv3_",k,"[AMRE_spl1$event==1]"))))
	tabRE <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMDG_spl1$","hv1_",k,"[AMDG_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMDG_spl1$","hv2_",k,"[AMDG_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMDG_spl1$","hv3_",k,"[AMDG_spl1$event==1]"))))
	tabDG <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMEC_spl1$","hv1_",k,"[AMEC_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMEC_spl1$","hv2_",k,"[AMEC_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMEC_spl1$","hv3_",k,"[AMEC_spl1$event==1]"))))
	tabEC <- rbind(tab1,tab2,tab3)
	
	tab1 <- table(eval(parse(text=paste0("AMOD_spl1$","hv1_",k,"[AMOD_spl1$event==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMOD_spl1$","hv2_",k,"[AMOD_spl1$event==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMOD_spl1$","hv3_",k,"[AMOD_spl1$event==1]"))))
	tabOD <- rbind(tab1,tab2,tab3)
	
	# Model for overall mortality - time-dependent #
	modTIb <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AM_spl"), 
		error=function(e) {NULL})

	AMH_spl1 <- AM_spl1[AM_spl1$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
	AMH_spl2 <- AM_spl2[AM_spl2$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
	AMH_spl3 <- AM_spl3[AM_spl3$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
	AMH_spl4 <- AM_spl4[AM_spl4$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
	AMH_spl5 <- AM_spl5[AM_spl5$charlsonSR==0,c("age","age_exit","out",paste0("hv1_",k),paste0("hv2_",k),paste0("hv3_",k))]
	
	tab1 <- table(eval(parse(text=paste0("AMH_spl1$","hv1_",k,"[AMH_spl1$out==1]"))))
	tab2 <- table(eval(parse(text=paste0("AMH_spl1$","hv2_",k,"[AMH_spl1$out==1]"))))
	tab3 <- table(eval(parse(text=paste0("AMH_spl1$","hv3_",k,"[AMH_spl1$out==1]"))))
	tabbH <- rbind(tab1,tab2,tab3)
	
		
	modTIbH <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMH_spl"), 
		error=function(e) {NULL})

	# Models for cause-specific mortality - time-dependent #
	modTICA <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMCA_spl"), 
		error=function(e) {NULL})

	modTICVD <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMCVD_spl"), 
		error=function(e) {NULL})

	modTIRE <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMRE_spl"), 
		error=function(e) {NULL})

	modTIDG <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMDG_spl"), 
		error=function(e) {NULL})

	modTIEC <- tryCatch(
coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMEC_spl"), 
		error=function(e) {NULL})

	modTIOD <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,event)~hv1_",k,"+hv2_",k,"+hv3_",k)), 5, N_MISSM, "AMOD_spl"), 
		error=function(e) {NULL})	
	
	# Measures to export #
		
	if(!is.null(modTIb))
	{
		FFb <- paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <-  tryCatch(pool(modTIb), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESb <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		} 				
	}	

	if(!is.null(modTIbH))
	{
		FFb <- paste0("Surv(age,age_exit,out)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <-  tryCatch(pool(modTIbH), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESbH <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		} 				
	}
		
	if(!is.null(modTICA))
	{
		FFCA <- paste0("Surv(age,age_exit,status==1)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTICA), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESCA <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}				
	}

	if(!is.null(modTICVD))
	{
		FFCVD <- paste0("Surv(age,age_exit,status==2)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTICVD), error=function(e) {NULL})	
		if (!is.null(pp))
		{
			RESCVD <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]), c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}				
	}

	if(!is.null(modTIRE))
	{
		FFRE <- paste0("Surv(age,age_exit,status==3)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIRE), error=function(e) {NULL})	
		if (!is.null(pp))
		{
			RESRE <- cbind(rep(levels(AM_spl1[,k]),3), c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}					
	}

	if(!is.null(modTIDG))
	{
		FFDG <- paste0("Surv(age,age_exit,status==4)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIDG), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESDG <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}			
	}

	if(!is.null(modTIEC))
	{
		FFEC <- paste0("Surv(age,age_exit,status==5)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIEC), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESEC <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]),c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))
		}		
	}

	if(!is.null(modTIOD))
	{
		FFOD <- paste0("Surv(age,age_exit,status==6)~hv1_",k,"+hv2_",k,"+hv3_",k)
		pp <- tryCatch(pool(modTIOD), error=function(e) {NULL})
		if (!is.null(pp))
		{
			RESOD <- cbind(rep(levels(AM_spl1[,k]),3),c("Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv1"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv2"],"Reference",pp$qbar[substr(names(pp$qbar),0,3)=="hv3"]), c("Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv1"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv2"],"Reference",sqrt(diag(pp$t))[substr(names(pp$qbar),0,3)=="hv3"]))	
		}					
	}	

	print(proc.time() - ptm)
	
	
	# Formulas
	YFFb <- rbind(YFFb,FFb)
	YFFCA <- rbind(YFFCA,FFCA)
	YFFCVD <- rbind(YFFCVD,FFCVD)
	YFFRE <- rbind(YFFRE,FFRE)
	YFFDG <- rbind(YFFDG,FFDG)
	YFFEC <- rbind(YFFEC,FFEC)
	YFFOD <- rbind(YFFOD,FFOD)

	# Coefficients + standard errors
	YRESb[[length(YRESb)+1]] <- RESb
	YRESCA[[length(YRESCA)+1]] <- RESCA
	YRESCVD[[length(YRESCVD)+1]] <- RESCVD
	YRESRE[[length(YRESRE)+1]] <- RESRE
	YRESDG[[length(YRESDG)+1]] <- RESDG
	YRESEC[[length(YRESEC)+1]] <- RESEC
	YRESOD[[length(YRESOD)+1]] <- RESOD

	YRESbH[[length(YRESbH)+1]] <- RESbH

	
	YTABb[[length(YTABb)+1]] <- tabb
	YTABCA[[length(YTABCA)+1]] <- tabCA
	YTABCVD[[length(YTABCVD)+1]] <- tabCVD
	YTABRE[[length(YTABRE)+1]] <- tabRE
	YTABDG[[length(YTABDG)+1]] <- tabDG
	YTABEC[[length(YTABEC)+1]] <- tabEC
	YTABOD[[length(YTABOD)+1]] <- tabOD
	
	YTABbH[[length(YTABbH)+1]] <- tabbH
	
	
	print(k)
}


## SAVE OBJECTS ##
for (i in c("YFFb","YFFCA","YFFCVD","YFFRE","YFFDG","YFFEC","YFFOD",
"YRESb","YRESCA","YRESCVD","YRESRE","YRESDG","YRESEC","YRESOD","YRESbH",
"YTABb","YTABCA","YTABCVD","YTABRE","YTABDG","YTABEC","YTABOD","YTABbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))}