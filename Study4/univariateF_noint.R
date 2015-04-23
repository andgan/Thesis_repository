#----------------------------------------------
# Filename: UnivariteF_noint.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
# Purpose: Univariate analysis in females. It obtaines both association results and C-index results for non-stratified analysis.
# Note: Also obtains the R^2 for association with age and the Schoenfeld's residuals P-values
#-----------------------------------------------
# Data used:  bdE9F, AF1-AF5.Rdata
# Data created: C-index: XCTbNI.Rdata,XCTCANI.Rdata,XCTCVDNI.Rdata,XCTRENI.Rdata,XCTDGNI.Rdata,XCTECNI.Rdata,XCTODNI.Rdata,XCTbHNI.Rdata
#								Formulas: XFFbNI.Rdata,XFFCANI.Rdata,XFFCVDNI.Rdata,XFFRENI.Rdata,XFFDGNI.Rdata,XFFECNI.Rdata,XFFODNI.Rdata,
#								Association results: XRESbNI.Rdata,XRESCANI.Rdata,XRESCVDNI.Rdata,XRESRENI.Rdata,XRESDGNI.Rdata,XRESECNI.Rdata,XRESODNI.Rdata,XRESbHNI.Rdata
#								Number of deaths per category: XTABbNI.Rdata,XTABCANI.Rdata,XTABCVDNI.Rdata,XTABRENI.Rdata,XTABDGNI.Rdata,XTABECNI.Rdata,XTABODNI.Rdata,XTABbHNI.Rdata
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

load("/proj/b2011036/uk.biobank/imputation_results/bdE9F.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

N_MISSF <- bdE9F$nmis


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

	# Calculate age - variable R2, this is useful for the final graph #
	r2 <- cor(AF1$age,predict(glm(AF1$age~AF1[,k])))^2
		
	# Determine the P-value for time-dependent effect #
	PROPp <- tryCatch({ 
	tb <- cox.zph(coxph(as.formula(paste0("Surv(age,age_exit,out==1)~",k)), data=AF1))
	f <- tb$table[nrow(tb$table),ncol(tb$table)]},
	error=function(e) {return(f)})
	PROPp <- ifelse(is.null(PROPp) | is.na(PROPp),1,PROPp)
		

	# Set reference level to most frequent level
	tt <- table(AF1[,k])
	AF1[,k]  <- relevel(AF1[,k], ref=names(tt)[tt==max(tt)][1])
	AF2[,k]  <- relevel(AF2[,k], ref=names(tt)[tt==max(tt)][1])
	AF3[,k]  <- relevel(AF3[,k], ref=names(tt)[tt==max(tt)][1])
	AF4[,k]  <- relevel(AF4[,k], ref=names(tt)[tt==max(tt)][1])
	AF5[,k]  <- relevel(AF5[,k], ref=names(tt)[tt==max(tt)][1])

	# N.of events in each class #
	tabb <- table(eval(parse(text=paste0("AF1$",k,"[AF1$out==1]"))))
	tabCA <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==1]"))))
	tabCVD <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==2]"))))
	tabRE <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==3]"))))
	tabDG <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==4]"))))
	tabEC <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==5]"))))
	tabOD <- table(eval(parse(text=paste0("AF1$",k,"[AF1$status==6]"))))
	
	
	# Model for overall mortality #
	modb <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~",k)), 5, N_MISSF,"AF"),
		error=function(e) {NULL})
		
	AFH1 <- AF1[AF1$charlsonSR==0,c("age","age_exit","out",k)]
	AFH2 <- AF2[AF2$charlsonSR==0,c("age","age_exit","out",k)]
	AFH3 <- AF3[AF3$charlsonSR==0,c("age","age_exit","out",k)]
	AFH4 <- AF4[AF4$charlsonSR==0,c("age","age_exit","out",k)]
	AFH5 <- AF5[AF5$charlsonSR==0,c("age","age_exit","out",k)]

	tabbH <- table(eval(parse(text=paste0("AFH1$",k,"[AFH1$out==1]"))))

	modbH <- tryCatch(
		coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~",k)), 5, N_MISSF,"AFH"),
		error=function(e) {NULL})	
					
	# Models for cause-specific mortality #
	modCA <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==1)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})

	modCVD <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==2)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})

	modRE <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==3)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})

	modDG <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==4)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})

	modEC <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==5)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})

	modOD <- tryCatch(
		coxph.mids(as.formula(paste("Surv(age,age_exit,status==6)~",k,sep="")), 5, N_MISSF,"AF"),
		error=function(e) {NULL})
	

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
						tryCatch(coxCC.temp1 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AFT1[-indCV,]),
						error=function(e) {NULL})
						tryCatch(coxCC.temp2 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AFT2[-indCV,]),
						error=function(e) {NULL})
						tryCatch(coxCC.temp3 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AFT3[-indCV,]),
						error=function(e) {NULL})
						tryCatch(coxCC.temp4 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AFT4[-indCV,]),
						error=function(e) {NULL})
						tryCatch(coxCC.temp5 <- CSC(as.formula(paste0("Hist(surv,status)~age+",k)),data=AFT5[-indCV,]),
						error=function(e) {NULL})
						
						CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"[indCV,]$surv"))), eval(parse(text=paste0("AFT",x,"[indCV,]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AFT",x,"[-indCV,]")))), newdata=eval(parse(text=paste0("AFT",x,"[indCV,]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))

						CTbH <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]$surv"))), eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AFT",x,"[-intersect(indCV,which(AFT1$charlsonSR==0)),]")))), newdata=eval(parse(text=paste0("AFT",x,"[intersect(indCV,which(AFT1$charlsonSR==0)),]")))))$concordance},mc.cores=5,mc.preschedule=FALSE))
											
											
						if(!is.null(coxCC.temp1))
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
						return(list(CTb,CTCA,CTCVD,CTRE,CTDG,CTEC,CTOD,CTbH))
						}



			
			
	# Calculate the measures to export
	
	if(!is.null(modb))
	{
			
		FFb <- paste0("Surv(age,age_exit,out)~",k)
		pp <- pool(modb)
		RESb <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))				
	}	
	
	if(!is.null(modbH))
	{
		pp <- pool(modbH)
		RESbH <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))
	}
	
	if(!is.null(modCA))
	{
		FFCA <- paste0("Surv(age,age_exit,status==1)~",k)
		pp <- pool(modCA)
		RESCA <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))			
	}
	
	if(!is.null(modCVD))
	{
		FFCVD <- paste0("Surv(age,age_exit,status==2)~",k)
		pp <- pool(modCVD)
		RESCVD <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))							
	}
	
	if(!is.null(modRE))
	{
		FFRE <- paste0("Surv(age,age_exit,status==3)~",k)
		pp <- pool(modRE)
		RESRE <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))					
	}
	
	if(!is.null(modDG))
	{
		FFDG <- paste0("Surv(age,age_exit,status==4)~",k)
		pp <- pool(modDG)
		RESDG <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))						
	}

	if(!is.null(modEC))
	{
		FFEC <- paste0("Surv(age,age_exit,status==5)~",k)
		pp <- pool(modEC)
		RESEC <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))					
	}
	
	if(!is.null(modOD))
	{
		FFOD <- paste0("Surv(age,age_exit,status==6)~",k)
		pp <- pool(modOD)
		RESOD <- cbind(levels(AF1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))				
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
	
	R2F <- c(R2F,r2)
	PROPpF <- c(PROPpF,PROPp)
	
	print(k)
}

## SAVE OBJECTS ##
for (i in c(
"XCTb","XCTCA","XCTCVD","XCTRE","XCTDG","XCTEC","XCTOD","XCTbH",
"XFFb","XFFCA","XFFCVD","XFFRE","XFFDG","XFFEC","XFFOD",
"XRESb","XRESCA","XRESCVD","XRESRE","XRESDG","XRESEC","XRESOD","XRESbH","R2F","PROPpF",
"XTABb","XTABCA","XTABCVD","XTABRE","XTABDG","XTABEC","XTABOD","XTABbH"))
{save(list=i, file=paste0("/proj/b2011036/uk.biobank/Univariate_results/no_int/",i,"NI.Rdata"))}