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
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

AF1 <- AF1[AF1$age>=40 & AF1$age<=70,]
AF2 <- AF2[AF2$age>=40 & AF2$age<=70,]
AF3 <- AF3[AF3$age>=40 & AF3$age<=70,]
AF4 <- AF4[AF4$age>=40 & AF4$age<=70,]
AF5 <- AF5[AF5$age>=40 & AF5$age<=70,]


## CREATE VALIDATION AND TRAINING DATASETS ##
## EXTERNAL VALIDATION: Glasgow and Edinburgh) ##
set.seed(123)
idsplitF2 <- which(AF1$center%in%c("Glasgow","Edinburgh"))

AFT1 <- AF1[-unique(c(idsplitF2)),]
AFV1 <- AF1[idsplitF2,]
AFT2 <- AF2[-unique(c(idsplitF2)),]
AFV2 <- AF2[idsplitF2,]
AFT3 <- AF3[-unique(c(idsplitF2)),]
AFV3 <- AF3[idsplitF2,]
AFT4 <- AF4[-unique(c(idsplitF2)),]
AFV4 <- AF4[idsplitF2,]
AFT5 <- AF5[-unique(c(idsplitF2)),]
AFV5 <- AF5[idsplitF2,]


NAMESF <- colnames(AF1)[!colnames(AF1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

resF <- read.xls("/proj/b2011036/uk.biobank/univariateF.xlsx", header=T, stringsAsFactor=F)

## BUILD THE FORMULA ACCORDING TO THE C-INDEX P-value ##
FFF <- NULL
for (f in NAMESF)
{		
	if (!any(resF$C.index.all[resF$Original.code == f] == "Singular fit") & !any(resF$C.index.all[resF$Original.code == f] == "N <= 5"))
	{
			k <- ifelse(unique(resF$P.value.for.Schoenfeld.residuals[resF$Original.code == f]) == "&lt;0.00001" , paste0(unique(resF$Original.code[resF$Original.code == f]),"*age"), unique(resF$Original.code[resF$Original.code == f]))
	}
	FFF <- c(FFF,k)
}


formufinF <- paste0("~",paste0(FFF, collapse="+"))


registerDoMC(cores=10)
for(k in 1:5)
{
	print(k)
	X <- model.matrix(as.formula(formufinF), eval(parse(text=paste0("AFT",k))))[,-1]

	set.seed(123)
	# Run penalized cox
	p.fac = rep(1, ncol(X))
	p.fac[which(colnames(X)=="age")] = 0
	
	lassomodFCV <- cv.glmnet(x=X, y=Surv(eval(parse(text=paste0("AFT",k,"$surv"))),eval(parse(text=paste0("AFT",k,"$out")))),family="cox", standardize=F,alpha=1, type.measure="deviance", nfolds=10, penalty.factor = p.fac, nlambda=100, parallel=TRUE)

	# find lambda for which mean squared error is minumum and the number of variables is between 10 and 20.
	ind  <- which(lassomodFCV$glmnet.fit$df < 8 & lassomodFCV$glmnet.fit$df> 3)
	optimal.lambda.index <- which(lassomodFCV$cvm==min(lassomodFCV$cvm[ind]))
	optimal.lambda  <- lassomodFCV$lambda[optimal.lambda.index]
	
	# find out the number of selected variables
	optimal.beta  <- lassomodFCV$glmnet.fit$beta[,optimal.lambda.index] 
	selectedBeta <- optimal.beta[abs(optimal.beta)>0 ] 
	
	# Find out the number of variables included in the true best model
	optimal.betaT  <- lassomodFCV$glmnet.fit$beta[,lassomodFCV$lambda==lassomodFCV$lambda.min] 
	selectedBetaT <- optimal.betaT[abs(optimal.betaT)>0 ] 
	
	# Calculate the C-index using Scotland
	XVAL <- model.matrix(as.formula(formufinF), data=eval(parse(text=paste0("AFV",k))))
	XVAL <- XVAL[,colnames(XVAL)%in%rownames(coef(lassomodFCV,s=optimal.lambda))]
	CCF <- survConcordance(Surv(eval(parse(text=paste0("AFV",k,"$surv"))),eval(parse(text=paste0("AFV",k,"$out"))))~predict(lassomodFCV,newx=XVAL,s=optimal.lambda))$concordance
	
	CCFT <- survConcordance(Surv(eval(parse(text=paste0("AFV",k,"$surv"))),eval(parse(text=paste0("AFV",k,"$out"))))~predict(lassomodFCV,newx=XVAL,s="lambda.min"))$concordance
	
	NNF <- length(selectedBeta)
	NNFT <- length(selectedBetaT)
	
	save(NNF,file=paste0("/proj/b2011036/uk.biobank/prediction_results/NNF",k,".Rdata"))
	save(NNFT,file=paste0("/proj/b2011036/uk.biobank/prediction_results/NNFT",k,".Rdata"))
	save(CCF,file=paste0("/proj/b2011036/uk.biobank/prediction_results/CCF",k,".Rdata"))
	save(CCFT,file=paste0("/proj/b2011036/uk.biobank/prediction_results/CCFT",k,".Rdata"))
}


