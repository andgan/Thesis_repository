#----------------------------------------------
# Filename: UnivariteM_noint.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in males. It obtaines both association results and C-index results for non-stratified analysis.
# Note: Also obtains the R^2 for association with age and the Schoenfeld's residuals P-values
#-----------------------------------------------
# Data used:  bdE9M.Rdata, AM1-AM5.Rdata
# Data created: C-index: YCTbNI.Rdata,YCTCANI.Rdata,YCTCVDNI.Rdata,YCTRENI.Rdata,YCTDGNI.Rdata,YCTECNI.Rdata,YCTODNI.Rdata,YCTbHNI.Rdata
#								Formulas: YFFbNI.Rdata,YFFCANI.Rdata,YFFCVDNI.Rdata,YFFRENI.Rdata,YFFDGNI.Rdata,YFFECNI.Rdata,YFFODNI.Rdata,
#								Association results: YRESbNI.Rdata,YRESCANI.Rdata,YRESCVDNI.Rdata,YRESRENI.Rdata,YRESDGNI.Rdata,YRESECNI.Rdata,YRESODNI.Rdata,YRESbHNI.Rdata
#								Number of deaths per category: YTABbNI.Rdata,YTABCANI.Rdata,YTABCVDNI.Rdata,YTABRENI.Rdata,YTABDGNI.Rdata,YTABECNI.Rdata,YTABODNI.Rdata,YTABbHNI.Rdata
#								Extra data: R2MNI.Rdata,PROPpMNI.Rdata

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

trsh <- 0.00001

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
	
	# Calculate age - variable R2, this is useful for the final graph #
	r2 <- cor(AM1$age,predict(glm(AM1$age~AM1[,k])))^2
	
	# Determine the P-value for time-dependent effect #
	PROPp <- tryCatch({ 
	tb <- cox.zph(coxph(as.formula(paste0("Surv(age,age_exit,out==1)~",k)), data=AM1))
	f <- tb$table[nrow(tb$table),ncol(tb$table)]},
	error=function(e) {return(f)})
	PROPp <- ifelse(is.null(PROPp) | is.na(PROPp),1,PROPp)
		

	# Set reference level to most frequent level
	tt <- table(AM1[,k])
	AM1[,k]  <- relevel(AM1[,k], ref=names(tt)[tt==max(tt)][1])
	AM2[,k]  <- relevel(AM2[,k], ref=names(tt)[tt==max(tt)][1])
	AM3[,k]  <- relevel(AM3[,k], ref=names(tt)[tt==max(tt)][1])
	AM4[,k]  <- relevel(AM4[,k], ref=names(tt)[tt==max(tt)][1])
	AM5[,k]  <- relevel(AM5[,k], ref=names(tt)[tt==max(tt)][1])

	# N.of events in each class #
	tabb <- table(eval(parse(text=paste0("AM1$",k,"[AM1$out==1]"))))
	tabCA <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==1]"))))
	tabCVD <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==2]"))))
	tabRE <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==3]"))))
	tabDG <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==4]"))))
	tabEC <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==5]"))))
	tabOD <- table(eval(parse(text=paste0("AM1$",k,"[AM1$status==6]"))))
	
	
	# Model for overall mortality #
	modb <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~",k)), 5, N_MISSM,"AM"),
		error=function(e) {NULL})
	
	# Healthy individuals (charlson score = 0)
	AMH1 <- AM1[AM1$charlsonSR==0,c("age","age_exit","out",k)]
	AMH2 <- AM2[AM2$charlsonSR==0,c("age","age_exit","out",k)]
	AMH3 <- AM3[AM3$charlsonSR==0,c("age","age_exit","out",k)]
	AMH4 <- AM4[AM4$charlsonSR==0,c("age","age_exit","out",k)]
	AMH5 <- AM5[AM5$charlsonSR==0,c("age","age_exit","out",k)]
	
	tabbH <- table(eval(parse(text=paste0("AMH1$",k,"[AMH1$out==1]"))))
		
	modbH <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~",k)), 5, N_MISSM,"AMH"),
		error=function(e) {NULL})	
		
	# Models for cause-specific mortality #
	modCA <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==1)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})

	modCVD <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==2)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})

	modRE <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==3)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})

	modDG <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==4)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})

	modEC <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==5)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})

	modOD <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==6)~",k,sep="")), 5, N_MISSM,"AM"),
		error=function(e) {NULL})
	
	# 10-fold C-index cross validation
	TEMP <- mclapply(1:10, function(i)
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

			indCV <- seq(i,nrow(AMT1),10)

			
			# Training set
			tryCatch(coxCC.temp1 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AMT1[-indCV,]),
			error=function(e) {NULL})
			tryCatch(coxCC.temp2 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AMT2[-indCV,]),
			error=function(e) {NULL})
			tryCatch(coxCC.temp3 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AMT3[-indCV,]),
			error=function(e) {NULL})
			tryCatch(coxCC.temp4 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AMT4[-indCV,]),
			error=function(e) {NULL})
			tryCatch(coxCC.temp5 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AMT5[-indCV,]),
			error=function(e) {NULL})

			# C-index overall mortality
			CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"[indCV,]$surv"))), eval(parse(text=paste0("AMT",x,"[indCV,]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AMT",x,"[-indCV,]")))), newdata=eval(parse(text=paste0("AMT",x,"[indCV,]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))
			
			# C-index healthy individuals
			CTbH <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]$surv"))), eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AMT",x,"[-intersect(indCV,which(AMT1$charlsonSR==0)),]")))), newdata=eval(parse(text=paste0("AMT",x,"[intersect(indCV,which(AMT1$charlsonSR==0)),]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))
			
			# C-index for cause-specific mortality
			if(!is.null(coxCC.temp1))
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
			return(list(CTb,CTCA,CTCVD,CTRE,CTDG,CTEC,CTOD,CTbH))
		}, mc.cores=16)	
	
	
	# Calculate the measures to export
	
	if(!is.null(modb))
	{
			
		FFb <- paste0("Surv(age,age_exit,out)~",k)
		pp <- pool(modb)
		RESb <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}	
	
	if(!is.null(modbH))
	{
		pp <- pool(modbH)
		RESbH <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
	
	if(!is.null(modCA))
	{
		FFCA <- paste0("Surv(age,age_exit,status==1)~",k)
		pp <- pool(modCA)
		RESCA <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
	
	if(!is.null(modCVD))
	{
		FFCVD <- paste0("Surv(age,age_exit,status==2)~",k)
		pp <- pool(modCVD)
		RESCVD <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))	
	}
	
	if(!is.null(modRE))
	{
		FFRE <- paste0("Surv(age,age_exit,status==3)~",k)
		pp <- pool(modRE)
		RESRE <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
	
	if(!is.null(modDG))
	{
		FFDG <- paste0("Surv(age,age_exit,status==4)~",k)
		pp <- pool(modDG)
		RESDG <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))	
	}

	if(!is.null(modEC))
	{
		FFEC <- paste0("Surv(age,age_exit,status==5)~",k)
		pp <- pool(modEC)
		RESEC <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
	
	if(!is.null(modOD))
	{
		FFOD <- paste0("Surv(age,age_exit,status==6)~",k)
		pp <- pool(modOD)
		RESOD <- cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
		
		
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


	R2M <- c(R2M,r2)
	PROPpM <- c(PROPpM,PROPp)
	
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
for (i in c(
"YCTb","YCTCA","YCTCVD","YCTRE","YCTDG","YCTEC","YCTOD","YCTbH",
"YFFb","YFFCA","YFFCVD","YFFRE","YFFDG","YFFEC","YFFOD",
"YRESb","YRESCA","YRESCVD","YRESRE","YRESDG","YRESEC","YRESOD","YRESbH","R2M","PROPpM",
"YTABb","YTABCA","YTABCVD","YTABRE","YTABDG","YTABEC","YTABOD","YTABbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/no_int/",i,"NI.Rdata"))}