#----------------------------------------------
# Filename: UnivariteF_yesint2.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in males. It obtaines C-index results for age-stratified analysis.
# Note: 
#-----------------------------------------------
# Data used:  AF1-AF5.Rdata
# Data created: C-index: XCTbYI.Rdata,XCTCAYI.Rdata,XCTCVDYI.Rdata,XCTREYI.Rdata,XCTDGYI.Rdata,XCTECYI.Rdata,XCTODYI.Rdata,XCTbHYI.Rdata
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

load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")


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
XCTb <- NULL
XCTCA <- NULL
XCTCVD <- NULL
XCTRE <- NULL
XCTDG <- NULL
XCTEC <- NULL
XCTOD <- NULL
XAIC <- NULL
XCTbH <- NULL

##########################
#### 3b. RUN ANALYSIS ####
##########################

for (k in NAMESF)
{	

		registerDoMC()

		TEMP <- foreach(i=1:10) %dopar%
			{

				  CTb <- -999
					CTCA <- -999
					CTCVD <- -999
					CTRE <- -999
					CTDG <- -999
					CTEC <- -999
					CTOD <- -999
					AIC <- -999
					CTbH <- -999

					indCV <- seq(i,nrow(AFT1),10)

					# Training set
					tryCatch(coxCC.temp1 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AFT1[-indCV,]),
					error=function(e) {NULL})
					tryCatch(coxCC.temp2 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AFT2[-indCV,]),
					error=function(e) {NULL})
					tryCatch(coxCC.temp3 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AFT3[-indCV,]),
					error=function(e) {NULL})
					tryCatch(coxCC.temp4 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AFT4[-indCV,]),
					error=function(e) {NULL})
					tryCatch(coxCC.temp5 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AFT5[-indCV,]),
					error=function(e) {NULL})


					CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"[indCV,]$surv"))), eval(parse(text=paste0("AFT",x,"[indCV,]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AFT",x,"[-indCV,]")))), newdata=eval(parse(text=paste0("AFT",x,"[indCV,]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))

					CTbH <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]$surv"))), eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AFT",x,"[-intersect(indCV,which(AFT1$charlsonSR==0)),]")))), newdata=eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))

					if (exists("coxCC.temp1") & exists("coxCC.temp2") & exists("coxCC.temp3") & exists("coxCC.temp4") & exists("coxCC.temp5"))
					{
						if(!is.null(coxCC.temp1) & !is.null(coxCC.temp2) & !is.null(coxCC.temp3) & !is.null(coxCC.temp4) & !is.null(coxCC.temp5))
						{
							cumHazAllCause1 <- mclapply(seq(1:6), 		function(c) { fastpredPARsample(coxCC.temp1$models[[paste("Cause", c)]], newdata = AFT1[indCV,],predtimes=5,ntosample=70, times=coxCC.temp1$eventTimes, n.cores=1)}, mc.cores=6)

							cumHazAllCause2 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp2$models[[paste("Cause", c)]], newdata = AFT2[indCV,],predtimes=5,ntosample=70, times=coxCC.temp2$eventTimes, n.cores=1)}, mc.cores=6)

							cumHazAllCause3 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp3$models[[paste("Cause", c)]], newdata = AFT3[indCV,],predtimes=5,ntosample=70, times=coxCC.temp3$eventTimes, n.cores=1)}, mc.cores=6)

							cumHazAllCause4 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp4$models[[paste("Cause", c)]], newdata = AFT4[indCV,],predtimes=5,ntosample=70, times=coxCC.temp4$eventTimes, n.cores=1)}, mc.cores=6)

							cumHazAllCause5 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp5$models[[paste("Cause", c)]], newdata = AFT5[indCV,],predtimes=5,ntosample=70, times=coxCC.temp5$eventTimes, n.cores=1)}, mc.cores=6)

							CTCA <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=1,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTCVD <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
								cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=2,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTRE <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
								cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=3,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTDG <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=4,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTEC <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=5,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTOD <- do.call("cbind", mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AFT",x,"[indCV,]"))) ,cause=6,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))
						}	
					}	
					return(list(CTb,CTCA,CTCVD,CTRE,CTDG,CTEC,CTOD,CTbH))
				}

	##  Collect the results ##

	print(TEMP)
	# C-statistics
	XCTb <- rbind(XCTb,colMeans(do.call("rbind",lapply(TEMP,"[[",1))))
	XCTCA <- rbind(XCTCA,colMeans(do.call("rbind",lapply(TEMP,"[[",2))))
	XCTCVD <- rbind(XCTCVD,colMeans(do.call("rbind",lapply(TEMP,"[[",3))))
	XCTRE <- rbind(XCTRE,colMeans(do.call("rbind",lapply(TEMP,"[[",4))))
	XCTDG <- rbind(XCTDG,colMeans(do.call("rbind",lapply(TEMP,"[[",5))))
	XCTEC <- rbind(XCTEC,colMeans(do.call("rbind",lapply(TEMP,"[[",6))))
	XCTOD <- rbind(XCTOD,colMeans(do.call("rbind",lapply(TEMP,"[[",7))))
	
	XCTbH <- rbind(XCTbH,colMeans(do.call("rbind",lapply(TEMP,"[[",8))))
		
	print(k)
}

## SAVE OBJECTS ##
for (i in c(
"XCTb","XCTCA","XCTCVD","XCTRE","XCTDG","XCTEC","XCTOD","XCTbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))}