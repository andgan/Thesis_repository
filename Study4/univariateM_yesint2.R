#----------------------------------------------
# Filename: UnivariteM_yesint2.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in males. It obtaines C-index results for age-stratified analysis.
# Note: 
#-----------------------------------------------
# Data used:  AM1-AM5.Rdata
# Data created: C-index: YCTbYI.Rdata,YCTCAYI.Rdata,YCTCVDYI.Rdata,YCTREYI.Rdata,YCTDGYI.Rdata,YCTECYI.Rdata,YCTODYI.Rdata,YCTbHYI.Rdata
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

load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")

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
YCTb <- NULL
YCTCA <- NULL
YCTCVD <- NULL
YCTRE <- NULL
YCTDG <- NULL
YCTEC <- NULL
YCTOD <- NULL
YAIC <- NULL
YCTbH <- NULL



##########################
#### 2b. RUN ANALYSIS ####
##########################

for (k in NAMESM)
{	

	
	ptm <- proc.time()
	
	TEMP <- mclapply(1:10, function(i)
		{
				CTb <- -999
				CTCA <- -999
				CTCVD <- -999
				CTRE <- -999
				CTDG <- -999
				CTEC <- -999
				CTOD <- -999
				CTbH <- -999				
			
			
				indCV <- seq(i,nrow(AMT1),10)

				# Training set
				tryCatch(coxCC.temp1 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AMT1[-indCV,]),
				error=function(e) {NULL})
				tryCatch(coxCC.temp2 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AMT2[-indCV,]),
				error=function(e) {NULL})
				tryCatch(coxCC.temp3 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AMT3[-indCV,]),
				error=function(e) {NULL})
				tryCatch(coxCC.temp4 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AMT4[-indCV,]),
				error=function(e) {NULL})
				tryCatch(coxCC.temp5 <- CSC(as.formula(paste0("Hist(surv,status)~age*",k)),data=AMT5[-indCV,]),
				error=function(e) {NULL})
				
				
				CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"[indCV,]$surv"))), eval(parse(text=paste0("AMT",x,"[indCV,]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x,"[-indCV,]")))), newdata=eval(parse(text=paste0("AMT",x,"[indCV,]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))
				
				CTbH <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]$surv"))), eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x,"[-intersect(indCV,which(AMT1$charlsonSR==0)),]")))), newdata=eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))
					
				if (exists("coxCC.temp1") & exists("coxCC.temp2") & exists("coxCC.temp3") & exists("coxCC.temp4") & exists("coxCC.temp5"))
				{
						if(!is.null(coxCC.temp1) & !is.null(coxCC.temp2) & !is.null(coxCC.temp3) & !is.null(coxCC.temp4) & !is.null(coxCC.temp5))
						{
							cumHazAllCause1 <- mclapply(seq(1:6), 		function(c) { fastpredPARsample(coxCC.temp1$models[[paste("Cause", c)]], newdata = AMT1[indCV,],predtimes=5,ntosample=70, times=coxCC.temp1$eventTimes, n.cores=4)}, mc.cores=6)

							cumHazAllCause2 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp2$models[[paste("Cause", c)]], newdata = AMT2[indCV,],predtimes=5,ntosample=70, times=coxCC.temp2$eventTimes, n.cores=4)}, mc.cores=6)

							cumHazAllCause3 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp3$models[[paste("Cause", c)]], newdata = AMT3[indCV,],predtimes=5,ntosample=70, times=coxCC.temp3$eventTimes, n.cores=4)}, mc.cores=6)

							cumHazAllCause4 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp4$models[[paste("Cause", c)]], newdata = AMT4[indCV,],predtimes=5,ntosample=70, times=coxCC.temp4$eventTimes, n.cores=4)}, mc.cores=6)

							cumHazAllCause5 <- mclapply(seq(1:6), 		function(c) {fastpredPARsample(coxCC.temp5$models[[paste("Cause", c)]], newdata = AMT5[indCV,],predtimes=5,ntosample=70, times=coxCC.temp5$eventTimes, n.cores=4)}, mc.cores=6)

							CTCA <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=1,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTCVD <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
								cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=2,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTRE <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
								cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=3,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTDG <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=4,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTEC <- do.call("cbind",mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=5,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))

							CTOD <- do.call("cbind", mclapply(1:5, function(x) {fastCindexsample(eval(parse(text=paste0("coxCC.temp",x))), formula=as.formula(paste0("Hist(surv,status)~1")),
							cens.model="marginal", data=eval(parse(text=paste0("AMT",x,"[indCV,]"))) ,cause=6,eval.times=5, causepred=get(paste0("cumHazAllCause",x)))}, mc.cores=5))
						}	
				}					
						return(list(CTb,CTCA,CTCVD,CTRE,CTDG,CTEC,CTOD,CTbH))
			}, mc.cores=16)
				
	# Measures to export #
		
	print(proc.time() - ptm)
	
	print(TEMP)
	##  Collect the results ##
	
	# C-statistics
	YCTb <- rbind(YCTb,colMeans(do.call("rbind",lapply(TEMP,"[[",1))))
	YCTCA <- rbind(YCTCA,colMeans(do.call("rbind",lapply(TEMP,"[[",2))))
	YCTCVD <- rbind(YCTCVD,colMeans(do.call("rbind",lapply(TEMP,"[[",3))))
	YCTRE <- rbind(YCTRE,colMeans(do.call("rbind",lapply(TEMP,"[[",4))))
	YCTDG <- rbind(YCTDG,colMeans(do.call("rbind",lapply(TEMP,"[[",5))))
	YCTEC <- rbind(YCTEC,colMeans(do.call("rbind",lapply(TEMP,"[[",6))))
	YCTOD <- rbind(YCTOD,colMeans(do.call("rbind",lapply(TEMP,"[[",7))))
	
	YCTbH <- rbind(YCTbH,colMeans(do.call("rbind",lapply(TEMP,"[[",8))))
	
	print(k)
}


## SAVE OBJECTS ##
for (i in c("YCTb","YCTCA","YCTCVD","YCTRE","YCTDG","YCTEC","YCTOD","YCTbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))}