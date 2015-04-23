library(foreign)
library(survival)
library(rms)
library(mice)
library(glmnet)
library(lattice)
library(WGCNA)
library(Hmisc)
library(foreach)
library(doMC)
library(doParallel)
library(mstate)
library(PredictABEL)
library(ggplot2)
library(MASS)
library(pec)
library(riskRegression)
library(cmprsk)
library(ResourceSelection)
library(gdata)

## Load in-house functions
source("/proj/b2011036/uk.biobank/Pgm/new_function.R")


## LOAD IMPUTED DATASETS ##
load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")

AM1 <- AM1[AM1$age>=40 & AM1$age<=70,]
AM2 <- AM2[AM2$age>=40 & AM2$age<=70,]
AM3 <- AM3[AM3$age>=40 & AM3$age<=70,]
AM4 <- AM4[AM4$age>=40 & AM4$age<=70,]
AM5 <- AM5[AM5$age>=40 & AM5$age<=70,]


## CREATE VALIDATION AND TRAINING DATASETS ##
## EXTERNAL VALIDATION: Glasgow and Edinburgh) ##
set.seed(123)
idsplitM2 <- which(AM1$center%in%c("Glasgow","Edinburgh"))

AMT1 <- AM1[-unique(c(idsplitM2)),]
AMV1 <- AM1[idsplitM2,]
AMT2 <- AM2[-unique(c(idsplitM2)),]
AMV2 <- AM2[idsplitM2,]
AMT3 <- AM3[-unique(c(idsplitM2)),]
AMV3 <- AM3[idsplitM2,]
AMT4 <- AM4[-unique(c(idsplitM2)),]
AMV4 <- AM4[idsplitM2,]
AMT5 <- AM5[-unique(c(idsplitM2)),]
AMV5 <- AM5[idsplitM2,]


NAMESM <- colnames(AM1)[!colnames(AM1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

resM <- read.xls("/proj/b2011036/uk.biobank/univariateM.xlsx", header=T, stringsAsFactor=F)

## BUILD THE FORMULA ACCORDING TO THE C-INDEX P-value ##
FFM <- NULL
for (f in NAMESM)
{		
	if (!any(resM$C.index.all[resM$Original.code == f] == "Singular fit") & !any(resM$C.index.all[resM$Original.code == f] == "N <= 5"))
	{
			k <- ifelse(unique(resM$P.value.for.Schoenfeld.residuals[resM$Original.code == f]) == "&lt;0.00001" , paste0(unique(resM$Original.code[resM$Original.code == f]),"*age"), unique(resM$Original.code[resM$Original.code == f]))
	}
	FFM <- c(FFM,k)
}


formufinM <- paste0("~",paste0(FFM, collapse="+"))


registerDoMC(cores=10)
for(k in 1:5)
{
	print(k)
	X <- model.matrix(as.formula(formufinM), eval(parse(text=paste0("AMT",k))))[,-1]

	set.seed(123)
	# Run penalized cox
	p.fac = rep(1, ncol(X))
	p.fac[which(colnames(X)=="age")] = 0
	
	lassomodMCV <- cv.glmnet(x=X, y=Surv(eval(parse(text=paste0("AMT",k,"$surv"))),eval(parse(text=paste0("AMT",k,"$out")))),family="cox", standardize=F,alpha=1, type.measure="deviance", nfolds=10, penalty.factor = p.fac, nlambda=100, parallel=TRUE)

	# find lambda for which mean squared error is minumum and the number of variables is between 10 and 20.
	ind  <- which(lassomodMCV$glmnet.fit$df < 8 & lassomodMCV$glmnet.fit$df> 3)
	optimal.lambda.index <- which(lassomodMCV$cvm==min(lassomodMCV$cvm[ind]))
	optimal.lambda  <- lassomodMCV$lambda[optimal.lambda.index]
	
	# find out the number of selected variables
	optimal.beta  <- lassomodMCV$glmnet.fit$beta[,optimal.lambda.index] 
	selectedBeta <- optimal.beta[abs(optimal.beta)>0 ] 
	
	# Find out the number of variables included in the true best model
	optimal.betaT  <- lassomodMCV$glmnet.fit$beta[,lassomodMCV$lambda==lassomodMCV$lambda.min] 
	selectedBetaT <- optimal.betaT[abs(optimal.betaT)>0 ] 
	
	# Calculate the C-index using Scotland
	XVAL <- model.matrix(as.formula(formufinM), data=eval(parse(text=paste0("AMV",k))))
	XVAL <- XVAL[,colnames(XVAL)%in%rownames(coef(lassomodMCV,s=optimal.lambda))]
	CC <- survConcordance(Surv(eval(parse(text=paste0("AMV",k,"$surv"))),eval(parse(text=paste0("AMV",k,"$out"))))~predict(lassomodMCV,newx=XVAL,s=optimal.lambda))$concordance
	
	CCT <- survConcordance(Surv(eval(parse(text=paste0("AMV",k,"$surv"))),eval(parse(text=paste0("AMV",k,"$out"))))~predict(lassomodMCV,newx=XVAL,s="lambda.min"))$concordance
	
	NN <- length(selectedBeta)
	NNT <- length(selectedBetaT)
	
	save(NN,file=paste0("/proj/b2011036/uk.biobank/prediction_results/NN",k,".Rdata"))
	save(NNT,file=paste0("/proj/b2011036/uk.biobank/prediction_results/NNT",k,".Rdata"))
	save(CC,file=paste0("/proj/b2011036/uk.biobank/prediction_results/CC",k,".Rdata"))
	save(CCT,file=paste0("/proj/b2011036/uk.biobank/prediction_results/CCT",k,".Rdata"))
}


