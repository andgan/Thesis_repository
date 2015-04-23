#----------------------------------------------
# Filename: Ana.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 16OCT014
# Updated: 02DEC2014 - Updated with additional analysis requested by reviewers
#					 15JAN2015 - Added new averages using variables from the census
# Purpose: Main analysis program for the UK biobank mortality project
# Note: This is the main program, which runs some other R scripts during execution: 
#				1. new_function.R: to load all the needed customized functions
#				2. code_variables_cat.R: to filter and cetegorized data
#       3. impute_uk_bioM.R, impute_uk_bioF.R: scripts to perform imputation, they are parallelized using 10 cores
#				4. UnivariateM_*.R and UnivariateF_*.R: scripts to calculate C-index and association result for each risk factor, both age-stratified and not.
#-----------------------------------------------
# Data used: check summary
# Data created: check summary
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

##############################
############
# SUMMARY:
# 1. LOAD AND PREPARE DATA. 
#		Data used: ukb2626.tab ukb3969.tab; 
#		Data created: out1.Rdata, out1bis, out2.Rdata, out3.Rdata, out4.Rdata, out5.Rdata, out7.Rdata, out8.Rdata. 
#		Pgm called: code_variables_cat.R

# 2. IMPUTATION: 
#		Data used: out8.Rdata, coded_questions.xlsx.
#	  Data created: out8M.Rdata, out8F.Rdata, pred_selM.Rdata, pred_selF.Rdata, bdE9M.Rdata, bdE9F.Rdata, AM1-5.Rdata, AF1-5.Rdata
#		Pgm called: impute_uk_bioM.R, impute_uk_bioF.R

# 3. DESCRIPTIVE STATISTICS:
#		Data used: AM1.Rdata, AF1.Rdata

# 4. UNIVARIATE ANALYSIS:
#		Data used: coded_questions.xlsx, AM1-5.Rdata, AF1-5.Rdata, bdE9M.Rdata, bdE9F.Rdata
#   Data created: data created by the scripts above (see each script for data created), univariateM.tab, univariateF.tab. These last two files are then manually edited and became: univariateM.xlsx, univariateF.xlsx
#	  Pgm called: UnivariateM_noint.R, UnivariateM_yesint1.R, UnivariateM_yesint2.R, UnivariateF_noint.R, UnivariateF_yesint1.R, UnivariateF_yesint2.R.

# 5. PREDICTION:
#	   Data used: AM1-5.Rdata, AF1-5.Rdata, bdE9M.Rdata, bdE9F.Rdata, univariateM.xlsx, univariateF.xlsx, coded_questions.xlsx
# 	 Data created: MAINMODM.Rdata, loeEM.Rdata, hlexpEM.Rdata, MAINMODF.Rdata, loeEF.Rdata, hlexpEF.Rdata

# 6. RE-WEIGHT BY MORTALITY TABLES AND CENSUS INFORMATION
#	   Data used: AM1-5.Rdata, AF1-5.Rdata, bdE9M.Rdata, bdE9F.Rdata, coded_questions.xlsx, general_health_M.csv, general_health_F.csv, van_car_M.csv, van_car_F.csv, male_2009_11.csv, female_2009_11.csv, MAINMODM.Rdata, MAINMODF.Rdata
# 	 Data created: GHMnewmean.Rdata, GHFnewmean.Rdata, VCMnewmean.Rdata, VCFnewmean.Rdata, wwMage.Rdata, wwFage.Rdata, MM095.Rdata, FF095.Rdata

# 7. CALCULATE THE RE-CALIBRATED INDIVIDUAL RISK AND SAVE MAIN QUANTITIES FOR PYTHON SCRIPT
#		 Data used: AM1-5.Rdata, AF1-5.Rdata, GHMnewmean.Rdata, GHFnewmean.Rdata, VCMnewmean.Rdata, VCFnewmean.Rdata, MAINMODM.Rdata, MAINMODF.Rdata, wwMage.Rdata, wwFage.Rdata
#		Data created:  varcovM.csv, varcovF.csv, coefmeanM.csv, coefmeanF.csv, basehazM.csv, basehazF.csv, All_riskM.Rdata, All_riskF.Rdata, AMV1EX.Rdata, AFV1EX.Rdata
		
# 8. BOOTSTRAP CI FOR CERTAIN QUANTITIES
#		 Data used: AM1-5.Rdata, AF1-5.Rdata

# 9. EXTRA ANALYSES FOR REVIWERS

#################################
#################################
### 1. LOAD AND PREPARE DATA ####
#################################
#################################


############################
### LOAD DATA - ID 2626 ####
############################

nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/uk.biobank/original_data/ukb2626.tab", sep=""), intern=T))," ")[[1]][1]
bd <- read.table("/proj/b2011036/uk.biobank/original_data/ukb2626.tab", header=T, comment.char = "", nrow=as.numeric(nrows), sep="\t")

## Import phenotypes ##
ph <- read.csv("/proj/b2011036/uk.biobank/var_and_codes.csv", header=T, stringsAsFactor=F)

## This runs all the labeling ##
source("/proj/b2011036/uk.biobank/original_data/ukb2626.r")

## Save dataset ##
save(bd,file="/proj/b2011036/uk.biobank/pre_imputation/out1.Rdata")

# N. of variables #
ncol(bd)
# N. of rows #
nrow(bd)

## Type of variables ##
table(unlist(lapply(bd,class)))


######################################################
### LOAD DATA - ID 3669 - UPDATED MORTALITY DATES ####
######################################################

nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/uk.biobank/original_data/ukb3969.tab", sep=""), intern=T))," ")[[1]][1]
bd2 <- read.table("/proj/b2011036/uk.biobank/original_data/ukb3969.tab", header=T, comment.char = "", nrow=as.numeric(nrows), sep="\t")

## This runs all the labeling ##
source("/proj/b2011036/uk.biobank/original_data/ukb3969.r")

## Save dataset ##
save(bd2,file="/proj/b2011036/uk.biobank/pre_imputation/out1bis.Rdata")

# N. of variables #
ncol(bd2)
# N. of rows #
nrow(bd2)

## Type of variables ##
table(unlist(lapply(bd2,class)))



############################
### MERGE TWO DATASETS  ####
############################
load("/proj/b2011036/uk.biobank/pre_imputation/out1bis.Rdata")
load("/proj/b2011036/uk.biobank/pre_imputation/out1.Rdata")


# Exclude these variables because a newer version is available in ID 3024
bd <- bd[,!colnames(bd)%in% c("f.40000.0.0","f.40000.1.0","f.40001.0.0","f.40001.1.0", 
 "f.40002.0.0",  "f.40002.0.1",  "f.40002.0.2",  "f.40002.0.3",  "f.40002.0.4", 
 "f.40002.0.5",  "f.40002.0.6",  "f.40002.0.7",  "f.40002.0.8",  "f.40002.0.9",
 "f.40002.0.10", "f.40002.0.11", "f.40002.0.12", "f.40002.0.13", "f.40002.1.0", 
 "f.40002.1.1",  "f.40002.1.2",  "f.40002.1.3",  "f.40002.1.4",  "f.40002.1.5", 
 "f.40002.1.6",  "f.40002.1.7",  "f.40002.1.8",  "f.40002.1.9",  "f.40002.1.10",
 "f.40002.1.11", "f.40002.1.12", "f.40002.1.13", "f.40007.0.0",  "f.40007.1.0", 
 "f.40018.0.0",  "f.40018.1.0")]

## Exclude all variables from pilot study (the first variable is kept since it the individual ID)##
bd <- bd[,c(TRUE,(unlist(lapply(strsplit(names(bd),"\\."),"[",3))!="1")[-1])]
sum(c(TRUE,(unlist(lapply(strsplit(names(bd),"\\."),"[",3))=="1")[-1]))


bdE1 <- merge(bd,bd2, by="f.eid")
save(bdE1,file="/proj/b2011036/uk.biobank/pre_imputation/out2.Rdata")

#######################
### CODE VARIABLES ####
#######################

## Exclude all partecipants from pilot study ##
bdE2 <- bdE1[bdE1$f.53.0.0 > "2007-01-01",]
sum(bdE1$f.53.0.0 < "2007-01-01")

## Exclude all those variables with 100% missing ##
excl_miss <- colSums(is.na(bdE2))
bdE3 <- bdE2[,excl_miss!=nrow(bdE2)]
names(bdE2)[excl_miss==nrow(bdE2)]
save(bdE3,file="/proj/b2011036/uk.biobank/pre_imputation/out3.Rdata")

## Define a new dataset ##
bdE4 <- bdE3

## TRANSFORM 'DON'T KNOW', 'UNSURE', 'PREFER NOT TO ANSWER' IN MISSING ##
for (i in 1:ncol(bdE4))
{
	if (class(bdE4[,i])=="ordered" | class(bdE4[,i])=="factor")
	{
		levels(bdE4[,i])[levels(bdE4[,i])%in%c("Prefer not to answer","Don't know","Unsure","Do not know")] <- NA
	}
	print(i)
}

save(bdE4,file="/proj/b2011036/uk.biobank/pre_imputation/out4.Rdata")


## RUN VARIABLE CODING ##
source("/proj/b2011036/uk.biobank/Pgm/code_variables_cat.R")

load("/proj/b2011036/uk.biobank/out5.Rdata")

table(unlist(lapply(bdE5,class)))
colnames(bdE5)[unlist(lapply(bdE5,class))=="numeric"]

#####################
## QUALITY CONTROL ##
#####################

## EXCLUDE SUBJECTS WITH >80% MISSING VARIABLES ##
subj_miss <- rowSums(is.na(bdE5))
bdE6 <- bdE5[(subj_miss/ncol(bdE5)) < 0.8,]
sum((subj_miss/ncol(bdE5)) >= 0.8)


## CORRECT SOME ISSUES, 2 SUBJECTS DO HAVE CAUSE OF DEATH, BUT NOT DEATH DATE ##
# Assign death at the median death date #
bdE6$ddate <- as.Date(bdE6$ddate, format="%Y-%m-%d")
bdE6$ddate[bdE6$out==0 & !is.na(bdE6$cofd)] <- median(bdE6$ddate, na.rm=T)

## DEFINE OUTCOMES, SURVIVAL TIME, AGE AT ENTRY AND EXIT ##
bdE6$surv <- bdE6$ddate-bdE6$date_center

## SCOTTISH CENTERAL HAVE A FOLLOW-UP UP TO 31DEC2012, ALL THE REST TO 17FEB2014
bdE6$surv[is.na(bdE6$surv) & !bdE6$center%in%c("Glasgow","Edinburgh")] <- max(bdE6$ddate,na.rm=T)-bdE6$date_center[is.na(bdE6$surv) & !bdE6$center%in%c("Glasgow","Edinburgh")]
bdE6$surv[is.na(bdE6$surv) & bdE6$center%in%c("Glasgow","Edinburgh")] <- as.Date("2012-12-31", format="%Y-%m-%d")-bdE6$date_center[is.na(bdE6$surv) & bdE6$center%in%c("Glasgow","Edinburgh")]

bdE6$surv <- as.numeric(bdE6$surv/365.25)
bdE6$age_exit <- bdE6$age + bdE6$surv
#bdE6$out <- ifelse(is.na(bdE6$ddate),0,1)
# Set max survival to 4 years
bdE6$out5 <- ifelse(bdE6$surv>5,0,bdE6$out)
bdE6$surv5 <- ifelse(bdE6$surv>5,5,bdE6$surv)
# Exclude subjects with negative survival or survival=0
bdE6 <- bdE6[bdE6$surv>0,] 
# Exclude these dates variable to avoid problems with imputations
bdE6$ddate <- NULL
bdE6$date_center <- NULL
# Code variable for competing risk #
bdE6$status <- ifelse(bdE6$cofd=="Cancer",1,
ifelse(bdE6$cofd=="Cardiovascular",2,
ifelse(bdE6$cofd=="Respiratory",3,
ifelse(bdE6$cofd=="Digestive",4,
ifelse(bdE6$cofd=="External causes",5,6)))))
bdE6$status[is.na(bdE6$status)] <- 0
# If no cause of death then is alive
bdE6$cofd <- as.factor(ifelse(is.na(bdE6$cofd), "Alive", as.character(bdE6$cofd)))


## PLOT N. MISSING PER VARIABLE ##
var_miss <- colSums(is.na(bdE6))
pdf("var_miss.pdf", width=12)
tp <- var_miss/nrow(bdE6)
barplot(sort(tp*100),axes=T,names.arg=names(tp)[order(tp)], cex.names=0.3,las=2, ylim=c(0,100),ylab="% missing")
abline(h=80, col="red")
dev.off()
# Name of variables with > 80% 
ph[as.character(ph[,3])%in%unlist(lapply(strsplit(names(bdE6)[(var_miss/nrow(bdE6))>0.80],"\\."),"[",2)),c(2,3)]

## EXCLUDE VARIABLES WITH > 80% MISSING SUBJECTS ##
bdE7 <- bdE6[,(var_miss/nrow(bdE6)) < 0.8]
sum((var_miss/nrow(bdE6)) >= 0.8)

save(bdE7,file="/proj/b2011036/uk.biobank/pre_imputation/out7.Rdata")

bdE8 <- bdE7


## EXCLUDE VARIABLES WITH ONLY ONE CLASS ##
excl_mis <- colnames(bdE8)[apply(bdE8,2,function(x) {length(unique(x[!is.na(x)]))==1})]
excl_mis
bdE8 <- bdE8[,!colnames(bdE8)%in%excl_mis]

## EXCLUDE VARIABLES WITH N-1 CLASSES < 20 DEATHS, with N=total NUMBER OF CLASSES ##
excl_less10 <- NULL
for (i in 1:ncol(bdE8))
{
	if (class(bdE8[,i])=="factor" | class(bdE8[,i])=="ordered")
	{
		ta <- table(bdE8[bdE8$out==1,i]) 
		te <- ifelse(sum(ta < 20) == length(ta)-1,colnames(bdE8)[i],NA) 
		excl_less10 <- cbind(excl_less10,te)
		print(i)
	}
}
excl_less10 <- excl_less10[!is.na(excl_less10)]
ph[as.character(ph[,3])%in%unlist(lapply(strsplit(excl_less10,"\\."),"[",2)),c(2,3)]

bdE8 <- bdE8[,!colnames(bdE8)%in%excl_less10]
# Name of the variables

save(bdE8,file="/proj/b2011036/uk.biobank/pre_imputation/out8.Rdata")


#### CHECK THAT THE TOP VARIABLES HAVE LOW MISSING VALUES #####
var_miss <- colSums(is.na(bdE8))
var_miss[colnames(bdE8)%in%MF]/nrow(bdE8)



###################
###################
## 2. IMPUTATION ##
###################
###################

# We created the Nelson-Aalen estimator to include in each imputation - see PMID: 19452569 #
hz <- basehaz(coxph(Surv(age,age_exit,out) ~ 1, data = bdE8))
idx <- match(bdE8$age_exit, hz[, "time"])
bdE8$nal <- hz[idx, "hazard"]

## SPECIFY SEX-SPECIFIC DATA ##
bdE8M <- bdE8[bdE8$sex=="Male",]
bdE8F <- bdE8[bdE8$sex=="Female",]

## EXCLUDE VARIABLES BELONGING TO THE OTHER SEX - 4041 excluded because only women, but coded in men as well ##
to_excM <- apply(bdE8M,2,function(x) sum(is.na(x))/length(x)<0.9)
bdE8M <- bdE8M[,c(colnames(bdE8M)[to_excM])]
bdE8M$f.4041.0.0 <- NULL

to_excF <- apply(bdE8F,2,function(x) sum(is.na(x))/length(x)<0.9)
bdE8F <- bdE8F[,c(colnames(bdE8F)[to_excF])]


save(bdE8M,file="/proj/b2011036/uk.biobank/pre_imputation/out8M.Rdata")
save(bdE8F,file="/proj/b2011036/uk.biobank/pre_imputation/out8F.Rdata")


## PERFORM SELECTION OF PREDICTORS ##
pred_selM <- quickpred2(bdE8M, mincor=0.1, minpuc=0, maxpred=10, nThreads=10, exclude=c("age","surv","age_exit","surv5","cofd","out5","status","sex"), include=c("out","nal","f.2178.0.0","center"))

## Check number of predictors per variable ##
rowSums(pred_selM)
rownames(pred_selM)[rowSums(pred_selM)==0]

save(pred_selM,file="/proj/b2011036/uk.biobank/imputation_results/pred_selM.Rdata")

pred_selF <- quickpred2(bdE8F, mincor=0.1, minpuc=0, maxpred=10, nThreads=10, exclude=c("age","surv","age_exit","surv5","cofd","out5","status","sex"), include=c("out","nal","f.2178.0.0","center"))

## Check number of predictors per variable ##
rowSums(pred_selF)
rownames(pred_selF)[rowSums(pred_selF)==0]

save(pred_selF,file="/proj/b2011036/uk.biobank/imputation_results/pred_selF.Rdata")

## IMPUTE - THIS IS DONE PARALLELIZED TO 10 CORES and using 512GB of RAM  ##
source("impute_uk_bioM.R")
source("impute_uk_bioF.R")


load("/proj/b2011036/uk.biobank/imputation_results/bdE9M.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/bdE9F.Rdata")


## READ DATA AND RBIND DATASETS FOR EACH IMPUTATION ##
AM1 <- complete(bdE9M,1)
AM2 <- complete(bdE9M,2)
AM3 <- complete(bdE9M,3)
AM4 <- complete(bdE9M,4)
AM5 <- complete(bdE9M,5)

# Correct a small error, there are 2 cases that have cause of death, but not overall death
AM1$out[AM1$out==0 & AM1$status!=0] <- 1
AM2$out[AM2$out==0 & AM2$status!=0] <- 1
AM3$out[AM3$out==0 & AM3$status!=0] <- 1
AM4$out[AM4$out==0 & AM4$status!=0] <- 1
AM5$out[AM5$out==0 & AM5$status!=0] <- 1
AM1$out5[AM1$out==0 & AM1$status!=0] <- 1
AM2$out5[AM2$out==0 & AM2$status!=0] <- 1
AM3$out5[AM3$out==0 & AM3$status!=0] <- 1
AM4$out5[AM4$out==0 & AM4$status!=0] <- 1
AM5$out5[AM5$out==0 & AM5$status!=0] <- 1

save(AM1,file="/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
save(AM2,file="/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
save(AM3,file="/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
save(AM4,file="/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
save(AM5,file="/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")


load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")


## FEMALE ##
AF1 <- complete(bdE9F,1)
AF2 <- complete(bdE9F,2)
AF3 <- complete(bdE9F,3)
AF4 <- complete(bdE9F,4)
AF5 <- complete(bdE9F,5)

save(AF1,file="/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
save(AF2,file="/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
save(AF3,file="/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
save(AF4,file="/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
save(AF5,file="/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

# Create variable with total number of missing #
N_MISSM <- bdE9M$nmis
N_MISSF <- bdE9F$nmis

# Check there are no missing data left #
pM <- apply(AM1, 2, function(x) sum(is.na(x)))
names(pM)[pM>0]
pF <- apply(AF1, 2, function(x) sum(is.na(x)))
names(pF)[pF>0]


### CHECK IF IMPUTED DATA ARE COMPARABLE WITH ORIGINAL DATA ###
longM <- AM1
longF <- AF1

missM <- is.na(bdE9M$data)
missF <- is.na(bdE9F$data)


# Only those variables with actual missing
codq <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)

pdf("/proj/b2011036/uk.biobank/chek_impM.pdf")	
for (i in colnames(longM[colnames(longM)%in%names(bdE9M$nmis[bdE9M$nmis!=0])]))
{	
	nameph <- codq[strsplit(i,"\\.")[[1]][2]==as.character(codq$Code),6]
	longt <- cbind(longM, hc.na=missM[,colnames(missM)==i])
	colnames(longt)[length(longt)] <- c("hc.na")
	noimp <- longt[longt$hc.na==F,i]
	yesimp <- longt[longt$hc.na==T,i]
	par(mar=c(10,2,2,2),oma = c(0, 0, 3, 0))
	xx <- barplot(cbind(table(yesimp)/length(yesimp)*100,table(noimp)/length(noimp)*100),axes=TRUE, names.arg=c(levels(noimp),levels(noimp)), beside=T, main="<- IMPUTED | NON IMPUTED ->", las=2,cex.names=0.8)
	text(x = xx[,1], y = 10, label = round(table(yesimp)/length(yesimp)*100,2), pos = 3, cex = 0.8, col = "red")
	text(x = xx[,1], y = 20, label = as.numeric(table(yesimp)), pos = 3, cex = 0.8, col = "blue")
	text(x = xx[,2], y = 10, label = round(table(noimp)/length(noimp)*100,2), pos = 3, cex = 0.8, col = "red")
	text(x = xx[,2], y = 20, label = as.numeric(table(noimp)), pos = 3, cex = 0.8, col = "blue")
	mtext(paste(nameph,"-",i), outer = TRUE ,cex=0.8)
	print(i)
}	
dev.off()


pdf("/proj/b2011036/uk.biobank/chek_impF.pdf")	
for (i in colnames(longF[colnames(longF)%in%names(bdE9F$nmis[bdE9F$nmis!=0])]))
{
	nameph <- codq[strsplit(i,"\\.")[[1]][2]==as.character(codq$Code),6]
	longt <- cbind(longF, hc.na=missF[,colnames(missF)==i])
	colnames(longt)[length(longt)] <- c("hc.na")
	noimp <- longt[longt$hc.na==F,i]
	yesimp <- longt[longt$hc.na==T,i]
	par(mar=c(10,2,2,2),oma = c(0, 0, 3, 0))
	xx <- barplot(cbind(table(yesimp)/length(yesimp)*100,table(noimp)/length(noimp)*100),axes=TRUE, names.arg=c(levels(noimp),levels(noimp)), beside=T, main="<- IMPUTED | NON IMPUTED ->", las=2,cex.names=0.8)
	text(x = xx[,1], y = 10, label = round(table(yesimp)/length(yesimp)*100,2), pos = 3, cex = 0.8, col = "red")
	text(x = xx[,1], y = 20, label = as.numeric(table(yesimp)), pos = 3, cex = 0.8, col = "blue")
	text(x = xx[,2], y = 10, label = round(table(noimp)/length(noimp)*100,2), pos = 3, cex = 0.8, col = "red")
	text(x = xx[,2], y = 20, label = as.numeric(table(noimp)), pos = 3, cex = 0.8, col = "blue")
	mtext(paste(nameph,"-",i), outer = TRUE ,cex=0.8)
	print(i)
}	
dev.off()


###################################
###################################
#### 3. DESCRIPTIVE STATISTICS ####
###################################
###################################

load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
nrow(AM1)
nrow(AF1)
mean(AM1$age)
sd(AM1$age)
mean(AF1$age)
sd(AF1$age)
table(AM1$out)
table(AF1$out)

summary(AM1$surv)
max(AM1$surv)
summary(AF1$surv)
max(AF1$surv)
table(AM1$status)
table(AF1$status)

load("/proj/b2011036/uk.biobank/out5.Rdata")
bdE4$sex <- bdE4$f.31.0.0
t2 <- substr(bdE4$f.40001.0.0,0,3)

head(table(t2[bdE4$cofd=="Cancer" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="Cancer" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="Cancer" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="Cancer" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$cofd=="Cardiovascular" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="Cardiovascular" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="Cardiovascular" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="Cardiovascular" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$cofd=="Respiratory" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="Respiratory" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="Respiratory" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="Respiratory" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$cofd=="Digestive" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="Digestive" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="Digestive" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="Digestive" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$cofd=="External causes" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="External causes" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="External causes" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="External causes" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$cofd=="other diseases" & bdE4$sex=="Male"])[order(-table(t2[bdE4$cofd=="other diseases" & bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$cofd=="other diseases" & bdE4$sex=="Female"])[order(-table(t2[bdE4$cofd=="other diseases" & bdE4$sex=="Female" ]))],5)

head(table(t2[bdE4$sex=="Male"])[order(-table(t2[bdE4$sex=="Male" ]))],5)
head(table(t2[bdE4$sex=="Female"])[order(-table(t2[bdE4$sex=="Female" ]))],5)



##################################
##################################
##### 4. UNIVARIATE ANALYSIS #####
##################################
##################################


## RUN UNIVARIATE SCRIPT (in parallel, 16 cores, 512GB ram) ##
# Unstratified analysis males, association and C-index results
source("UnivariateM_noint.R")
# Age-stratified analysis males, association results
source("UnivariateM_yesint1.R")
# Age-stratified analysis males, C-index results
source("UnivariateM_yesint2.R")
# Unstratified analysis females, association and C-index results
source("UnivariateF_noint.R")
# Age-stratified analysis females, association results
source("UnivariateF_yesint1.R")
# Age-stratified analysis females, C-index results
source("UnivariateF_yesint2.R")

##########
## MALE ##
##########


## LOAD DATA -  STRATIFIED ##
for (i in c(
"YCTb","YCTCA","YCTCVD","YCTRE","YCTDG","YCTEC","YCTOD","YCTbH",
"YFFb","YFFCA","YFFCVD","YFFRE","YFFDG","YFFEC","YFFOD",
"YRESb","YRESCA","YRESCVD","YRESRE","YRESDG","YRESEC","YRESOD","YRESbH",
"YTABb","YTABCA","YTABCVD","YTABRE","YTABDG","YTABEC","YTABOD","YTABbH"))
{load(paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))
assign(paste0(i,"YI"),get(i))}



## LOAD DATA - NON STRATIFIED ##
for (i in c(
"YCTb","YCTCA","YCTCVD","YCTRE","YCTDG","YCTEC","YCTOD","YCTbH",
"YFFb","YFFCA","YFFCVD","YFFRE","YFFDG","YFFEC","YFFOD",
"YRESb","YRESCA","YRESCVD","YRESRE","YRESDG","YRESEC","YRESOD","R2M","PROPpM","YRESbH",
"YTABb","YTABCA","YTABCVD","YTABRE","YTABDG","YTABEC","YTABOD","YTABbH"))
{load(paste0("/proj/b2011036/uk.biobank/Univariate_results/no_int/",i,"NI.Rdata"))
assign(paste0(i,"NI"),get(i))}



NAMESM <- colnames(AM1)[!colnames(AM1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

cutpoint <- c(53,62)


## PROCESS OBJECTS TO CREATE A FINAL LIST IN EXCEL FORMAT ##
# Assuming that YRESCA is without -999, used as reference
FINM <- NULL
for (dis in c("b","bH","CA","CVD","RE","DG","EC","OD"))
{
	TM <- NULL
	for (i in 1:length(get(paste0("YRES",dis,"NI"))))
	{
		
		pp <- ifelse(PROPpMNI[i]< 0.00001,"<0.00001",round(PROPpMNI[i],5))
		## NON STRATIFIED ##
			t <- eval(parse(text=paste0("YRES",dis,"NI[[i]]")))
			if (t == -999)
			{
				if (length(grep("l",NAMESM[i]))>0)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else if (strsplit(NAMESM[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESM[i],"\\.")[[1]][2]}
				mm <- matrix(c(NA,NA,NA,NA,NA,NA,extN,NA,R2M[i],pp,NAMESM[i]), nrow=1)
				tMNI <- mm[rep(1:nrow(mm), times = nrow(YRESCANI[[i]])), ]
				print(dis)
				print(i)
			}
			else
			{
				# This is to calculate the number of events in each class
				tab1 <- eval(parse(text=paste0("YTAB",dis,"NI")))[[i]]
				agetres <- "All"

				# This is to extract the variable code for subsequent matching
				if (length(grep("l",NAMESM[i]))>0)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else if (strsplit(NAMESM[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESM[i],"\\.")[[1]][2]}

				# Calculate P-value
				pval <- 2*pnorm(-abs(as.numeric(t[,2])/as.numeric(t[,3])))
				
				
				# Average C-index
				cindAV <- round(mean(eval(parse(text=paste0("YCT",dis,"NI[i,]")))),7)

				# Final cbind
				tMNI <- cbind(t,tab1,pval,agetres,extN,cindAV,R2M[i],pp,NAMESM[i])
				if (length(tMNI[as.numeric(tMNI[,4]) <= 5,4])== (nrow(tMNI)-1))
				{
					tMNI[,8] <- "N <= 5"
				}			
			}
				
		## AGE-STRATIFIED ##
			t <- eval(parse(text=paste0("YRES",dis,"YI[[i]]")))	
			## IF MISSING VALUES ##
			if (t == -999)
			{
				if (length(grep("l",NAMESM[i]))>0)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else if (strsplit(NAMESM[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESM[i],"\\.")[[1]][2]}
				mm <- matrix(c(NA,NA,NA,NA,NA,NA,extN,NA,R2M[i],pp,NAMESM[i]), nrow=1)
				tMYI <- mm[rep(1:nrow(mm), times = nrow(YRESCANI[[i]])*3), ]
				print(dis)
				print(i)
			}
			else
			{
				# This is to calculate the number of events in each class
				tab1 <- eval(parse(text=paste0("YTAB",dis,"YI")))[[i]][1,]
				tab2 <- eval(parse(text=paste0("YTAB",dis,"YI")))[[i]][2,]
				tab3 <- eval(parse(text=paste0("YTAB",dis,"YI")))[[i]][3,]

				# Paste the age class
				agetres <- c(rep(paste0("<",cutpoint[1]),length(tab1)) , rep(paste0(cutpoint[1],"-",cutpoint[2]),length(tab2)), rep(paste0(">",cutpoint[2]),length(tab3)))

				# This is to extract the variable code for subsequent matching
				if (length(grep("l",NAMESM[i]))>0)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])}else if (strsplit(NAMESM[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESM[i],"\\.")[[1]][2],"_",strsplit(NAMESM[i],"\\.")[[1]][4])}else {extN <- strsplit(NAMESM[i],"\\.")[[1]][2]}

				# Calculate P-value
				pval <- 2*pnorm(-abs(as.numeric(t[,2])/as.numeric(t[,3])))

				# Average C-index
				cindAV <- round(mean(eval(parse(text=paste0("YCT",dis,"YI[i,]")))),7)
				# Final cbind
				tMYI <- cbind(t,c(tab1,tab2,tab3),pval,agetres,extN,cindAV,R2M[i],pp,NAMESM[i])
				if (length(tMYI[as.numeric(tMYI[,4]) <= 5,4])==nrow(tMYI)-3)
				{
					tMYI[,8] <- "N <= 5"
				}			
			}
		
		TM <- rbind(TM,rbind(tMNI,tMYI))
	}
	
	# Report Hazard ratios
	TM[TM[,2]!="Reference" & !is.na(TM[,2]),2] <- paste0(formatC(exp(as.numeric(TM[TM[,2]!="Reference" & !is.na(TM[,2]),2])),digits=1, format="f") ," [",formatC(exp(as.numeric(TM[TM[,2]!="Reference" & !is.na(TM[,2]),2]) - 1.96*as.numeric(TM[TM[,3]!="Reference" & !is.na(TM[,3]),3])),digits=1, format="f"),"-",formatC(exp(as.numeric(TM[TM[,2]!="Reference" & !is.na(TM[,2]),2]) + 1.96*as.numeric(TM[TM[,3]!="Reference" & !is.na(TM[,3]),3])),digits=1, format="f"),"]")
	# Standard errors away
	TM <- TM[,-3] 
	# If sample size is too small then report N < 5 #
	TM[as.numeric(TM[,3]) <= 5,2] <- "N <= 5"
	TM[as.numeric(TM[,3]) <= 5,4] <- "N <= 5"
	# P-values as reference
	TM[TM[,2]=="Reference",4] <- "Reference"
	TM[TM[,2]=="Reference",5] <- "Reference"
	
	
	FINM <- cbind(FINM,TM)
}

# Remove some columns because duplicated #
FINM <- cbind(FINM[,c(6,10,8,9,5,1,2,3,4,7,12,13,14,17,22,23,24,27,32,33,34,37,42,43,44,47,52,53,54,57,62,63,64,67,72,73,74,77)],1:nrow(FINM))
colnames(FINM) <- c("UK Biobank code","Original code","R2 for association with age","P-value for Schoenfeld residuals", "Age category","Categories","Hazard Ratio [95% C.I.] all","N.deaths all","P-value all",  "C-index all","Hazard Ratio [95% C.I.] Healthy", "N.deaths Healthy","P-value Healthy", "C-index Healthy", "Hazard Ratio [95% C.I.] CA", "N.deaths CA","P-value CA", "C-index CA", "Hazard Ratio [95% C.I.] CVD","N.deaths CVD", "P-value CVD", "C-index CVD", "Hazard Ratio [95% C.I.] RE","N.deaths RE", "P-value RE", "C-index RE", "Hazard Ratio [95% C.I.] DG", "N.deaths DG", "P-value DG", "C-index DG", "Hazard Ratio [95% C.I.] EC","N.deaths EC", "P-value EC", "C-index EC", "Hazard Ratio [95% C.I.] OD", "N.deaths OD", "P-value OD", "C-index OD","uniqueID")


# Now merge with the file that includes the questions #
codq <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)

# Check all codes are also in the FIN file #
setdiff(FINM[,1],codq$Code)

# Merge
to_expM <- merge(FINM, codq, by.x="UK Biobank code", by.y="Code")
to_expM2 <- to_expM[order(as.numeric(as.character(to_expM$uniqueID))),]


write.table(to_expM2,file="/proj/b2011036/uk.biobank/univariateM.tab", quote=F, row.names=F, sep="$")




############
## FEMALE ##
############


## LOAD DATA -  STRATIFIED ##
for (i in c(
"XCTb","XCTCA","XCTCVD","XCTRE","XCTDG","XCTEC","XCTOD","XCTbH",
"XFFb","XFFCA","XFFCVD","XFFRE","XFFDG","XFFEC","XFFOD",
"XRESb","XRESCA","XRESCVD","XRESRE","XRESDG","XRESEC","XRESOD","XRESbH",
"XTABb","XTABCA","XTABCVD","XTABRE","XTABDG","XTABEC","XTABOD","XTABbH"))
{load(paste0("/proj/b2011036/uk.biobank/Univariate_results/yes_int/",i,"YI.Rdata"))
assign(paste0(i,"YI"),get(i))}


## LOAD DATA - NON STRATIFIED ##
for (i in c(
"XCTb","XCTCA","XCTCVD","XCTRE","XCTDG","XCTEC","XCTOD","XCTbH",
"XFFb","XFFCA","XFFCVD","XFFRE","XFFDG","XFFEC","XFFOD",
"XRESb","XRESCA","XRESCVD","XRESRE","XRESDG","XRESEC","XRESOD","R2F","PROPpF","XRESbH",
"XTABb","XTABCA","XTABCVD","XTABRE","XTABDG","XTABEC","XTABOD","XTABbH"))
{load(paste0("/proj/b2011036/uk.biobank/Univariate_results/no_int/",i,"NI.Rdata"))
assign(paste0(i,"NI"),get(i))}


NAMESF <- colnames(AF1)[!colnames(AF1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

cutpoint <- c(53,62)


## PROCESS OBJECTS TO CREATE A FINAL LIST IN EXCEL FORMAT ##

FINF <- NULL
for (dis in c("b","bH","CA","CVD","RE","DG","EC","OD"))
{
	TF <- NULL
	for (i in 1:length(get(paste0("XRES",dis,"NI"))))
	{
		
		pp <- ifelse(PROPpFNI[i]< 0.00001,"<0.00001",round(PROPpFNI[i],5))
		## NON STRATIFIED ##
			t <- eval(parse(text=paste0("XRES",dis,"NI[[i]]")))
			if (t == -999)
			{
				if (length(grep("l",NAMESF[i]))>0)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else if (strsplit(NAMESF[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESF[i],"\\.")[[1]][2]}
				mm <- matrix(c(NA,NA,NA,NA,NA,NA,extN,NA,R2F[i],pp,NAMESF[i]), nrow=1)
				tFNI <- mm[rep(1:nrow(mm), times = nrow(XRESCANI[[i]])), ]
				print(dis)
				print(i)
			}
			else
			{
				# This is to calculate the number of events in each class
				tab1 <- eval(parse(text=paste0("XTAB",dis,"NI")))[[i]]
				agetres <- "All"

				# This is to extract the variable code for subsequent matching
				if (length(grep("l",NAMESF[i]))>0)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else if (strsplit(NAMESF[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESF[i],"\\.")[[1]][2]}

				# Calculate P-value
				pval <- 2*pnorm(-abs(as.numeric(t[,2])/as.numeric(t[,3])))
				
				
				# Average C-index
				cindAV <- round(mean(eval(parse(text=paste0("XCT",dis,"NI[i,]")))),7)

				# Final cbind
				tFNI <- cbind(t,tab1,pval,agetres,extN,cindAV,R2F[i],pp,NAMESF[i])	
				if (length(tFNI[as.numeric(tFNI[,4]) <= 5,4])== (nrow(tFNI)-1))
				{
					tFNI[,8] <- "N <= 5"
				}			
			}
				
		## AGE-STRATIFIED ##
			t <- eval(parse(text=paste0("XRES",dis,"YI[[i]]")))	
			## IF MISSING VALUES ##
			if (t == -999)
			{
				if (length(grep("l",NAMESF[i]))>0)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else if (strsplit(NAMESF[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])} else {extN <- strsplit(NAMESF[i],"\\.")[[1]][2]}
				mm <- matrix(c(NA,NA,NA,NA,NA,NA,extN,NA,R2F[i],pp,NAMESF[i]), nrow=1)
				tFYI <- mm[rep(1:nrow(mm), times = nrow(XRESCANI[[i]])*3), ]
				print(dis)
				print(i)
			}
			else
			{
				# This is to calculate the number of events in each class
				tab1 <- eval(parse(text=paste0("XTAB",dis,"YI")))[[i]][1,]
				tab2 <- eval(parse(text=paste0("XTAB",dis,"YI")))[[i]][2,]
				tab3 <- eval(parse(text=paste0("XTAB",dis,"YI")))[[i]][3,]

				# Paste the age class
				agetres <- c(rep(paste0("<",cutpoint[1]),length(tab1)) , rep(paste0(cutpoint[1],"-",cutpoint[2]),length(tab2)), rep(paste0(">",cutpoint[2]),length(tab3)))

				# This is to extract the variable code for subsequent matching
				if (length(grep("l",NAMESF[i]))>0)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])}else if (strsplit(NAMESF[i],"\\.")[[1]][4]==2)
				{extN <- paste0(strsplit(NAMESF[i],"\\.")[[1]][2],"_",strsplit(NAMESF[i],"\\.")[[1]][4])}else {extN <- strsplit(NAMESF[i],"\\.")[[1]][2]}

				# Calculate P-value
				pval <- 2*pnorm(-abs(as.numeric(t[,2])/as.numeric(t[,3])))

				# Average C-index
				cindAV <- round(mean(eval(parse(text=paste0("XCT",dis,"YI[i,]")))),7)
				# Final cbind
				tFYI <- cbind(t,c(tab1,tab2,tab3),pval,agetres,extN,cindAV,R2F[i],pp,NAMESF[i])					
				if (length(tFYI[as.numeric(tFYI[,4]) <= 5,4])==nrow(tFYI)-3)
				{
					tFYI[,8] <- "N <= 5"
				}
			}
		
		TF <- rbind(TF,rbind(tFNI,tFYI))
	}
	
	# Report Hazard ratios
	TF[TF[,2]!="Reference" & !is.na(TF[,2]),2] <- paste0(formatC(exp(as.numeric(TF[TF[,2]!="Reference" & !is.na(TF[,2]),2])),digits=1, format="f") ," [",formatC(exp(as.numeric(TF[TF[,2]!="Reference" & !is.na(TF[,2]),2]) - 1.96*as.numeric(TF[TF[,3]!="Reference" & !is.na(TF[,3]),3])),digits=1, format="f"),"-",formatC(exp(as.numeric(TF[TF[,2]!="Reference" & !is.na(TF[,2]),2]) + 1.96*as.numeric(TF[TF[,3]!="Reference" & !is.na(TF[,3]),3])),digits=1, format="f"),"]")
	# Standard errors away
	TF <- TF[,-3] 
	# If sample size is too small then report N < 5 #
	TF[as.numeric(TF[,3]) <= 5,2] <- "N <= 5"
	TF[as.numeric(TF[,3]) <= 5,4] <- "N <= 5"
	# P-values as reference
	TF[TF[,2]=="Reference",4] <- "Reference"
	TF[TF[,2]=="Reference",5] <- "Reference"
	
	FINF <- cbind(FINF,TF)
}

# Remove some columns because duplicated #
FINF <- cbind(FINF[,c(6,10,8,9,5,1,2,3,4,7,12,13,14,17,22,23,24,27,32,33,34,37,42,43,44,47,52,53,54,57,62,63,64,67,72,73,74,77)],1:nrow(FINF))
colnames(FINF) <- c("UK Biobank code","Original code","R2 for association with age","P-value for Schoenfeld residuals", "Age category","Categories","Hazard Ratio [95% C.I.] all","N.deaths all","P-value all",  "C-index all","Hazard Ratio [95% C.I.] Healthy", "N.deaths Healthy","P-value Healthy", "C-index Healthy", "Hazard Ratio [95% C.I.] CA", "N.deaths CA","P-value CA", "C-index CA", "Hazard Ratio [95% C.I.] CVD","N.deaths CVD", "P-value CVD", "C-index CVD", "Hazard Ratio [95% C.I.] RE","N.deaths RE", "P-value RE", "C-index RE", "Hazard Ratio [95% C.I.] DG", "N.deaths DG", "P-value DG", "C-index DG", "Hazard Ratio [95% C.I.] EC","N.deaths EC", "P-value EC", "C-index EC", "Hazard Ratio [95% C.I.] OD", "N.deaths OD", "P-value OD", "C-index OD","uniqueID")


# Now merge with the file that includes the questions #
codq <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)

# Check all codes are also in the FIN file #
setdiff(FINF[,1],codq$Code)

# Merge
to_expF <- merge(FINF, codq, by.x="UK Biobank code", by.y="Code")
to_expF2 <- to_expF[order(as.numeric(as.character(to_expF$uniqueID))),]


write.table(to_expF2,file="/proj/b2011036/uk.biobank/univariateF.tab", quote=F, row.names=F, sep="$")




#########################
#########################
##### 5. PREDICTION #####
#########################
#########################


## LOAD IMPUTED DATASETS ##
load("/proj/b2011036/uk.biobank/imputation_results/bdE9M.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AM5.Rdata")

N_MISSM <- bdE9M$nmis


load("/proj/b2011036/uk.biobank/imputation_results/bdE9F.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF2.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF3.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF4.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF5.Rdata")

N_MISSF <- bdE9F$nmis

## FOR THE PREDICTION FROM NOW-ON WE EXCLUDE THOSE <40 and >70 BECAUSE TO FEW CASES AND
## BECAUSE OUR PREDICTION DO NOT WANT TO EXAPAND OUT OF THAT RANGE ##

AM1 <- AM1[AM1$age>=40 & AM1$age<=70,]
AM2 <- AM2[AM2$age>=40 & AM2$age<=70,]
AM3 <- AM3[AM3$age>=40 & AM3$age<=70,]
AM4 <- AM4[AM4$age>=40 & AM4$age<=70,]
AM5 <- AM5[AM5$age>=40 & AM5$age<=70,]

AF1 <- AF1[AF1$age>=40 & AF1$age<=70,]
AF2 <- AF2[AF2$age>=40 & AF2$age<=70,]
AF3 <- AF3[AF3$age>=40 & AF3$age<=70,]
AF4 <- AF4[AF4$age>=40 & AF4$age<=70,]
AF5 <- AF5[AF5$age>=40 & AF5$age<=70,]

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




###########
### MALE ##
###########

NAMESM <- colnames(AM1)[!colnames(AM1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

resM <- read.xls("/proj/b2011036/uk.biobank/univariateM.xlsx", header=T, stringsAsFactor=F)

##EXCLUDE ALL CATEGORIES SHOULD NOT BE INCLUDED IN THE SCORE ##
resM2 <- resM[!resM$Measurement.class%in%c("Physical measures","Blood assays","Cognitive function"),]


## EXTRACT TOP 20 MEASUREMENTS FOR EACH CAUSE OF DEATH
FFM20 <- NULL
for (dis in c("CA","CVD","RE","DG","EC","OD"))
{
	K <- NULL
	for (f in unique(resM2$Original.code))
	{
		
		t <- resM2[resM2$Original.code == f,]
		if(unique(t$P.value.for.Schoenfeld.residuals)!= "&lt;0.00001")
		{
			t2 <- t[t$Age.category=="All",]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),unique(t$Original.code),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}	
					else {k <- c(unique(t$Original.code),unique(t$Original.code),NA)}
							
		}	
		else if(unique(t$P.value.for.Schoenfeld.residuals)== "&lt;0.00001")
		{
			t2 <- t[!t$Age.category%in%c("All","Reference"),]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}		
					else {k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),NA)}		
		}
		
		K <- rbind(K,k)
	}	
	
	FFM20 <- cbind(FFM20,K[order(-as.numeric(K[,3])),2][1:20])
}



to_sel <- unique(as.vector(FFM20))
formufinM <- paste0("~",paste0(to_sel, collapse="+"))



## FORMULA INCLUDING ALL THE FACTORS (THIS IS USED LATER) ##
FFM <- NULL
for (f in NAMESM)
{		
	if (!any(resM$C.index.all[resM$Original.code == f] == "Singular fit") & !any(resM$C.index.all[resM$Original.code == f] == "N <= 5"))
	{
			k <- ifelse(unique(resM$P.value.for.Schoenfeld.residuals[resM$Original.code == f]) == "&lt;0.00001" , paste0(unique(resM$Original.code[resM$Original.code == f]),"*age"), unique(resM$Original.code[resM$Original.code == f]))
	}
	FFM <- c(FFM,k)
}



###############################################################
## LASSO REGRESSION - RUN PARALLELIZED TO SPEED CALCULATIONS ##
###############################################################

# Formula for lasso, including all the measurements
formufinML <- paste0("~",paste0(FFM, collapse="+"))

registerDoMC(cores=10)
CCM <- NULL
NNM <- NULL
CCMT <- NULL
NNMT <- NULL	
for(k in 1:5)
{
	print(k)
	X <- model.matrix(as.formula(formufinML), eval(parse(text=paste0("AMT",k))))[,-1]

	set.seed(123)
	# Run penalized cox
	p.fac = rep(1, ncol(X))
	p.fac[which(colnames(X)=="age")] = 0
	
	lassomodMCV <- cv.glmnet(x=X, y=Surv(eval(parse(text=paste0("AMT",k,"$surv"))),eval(parse(text=paste0("AMT",k,"$out")))),family="cox", standardize=F,alpha=1, type.measure="deviance", nfolds=10, penalty.factor = p.fac, nlambda=200, parallel=TRUE)

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
	XVAL <- model.matrix(as.formula(formufinML), data=eval(parse(text=paste0("AMV",k))))
	XVAL <- XVAL[,colnames(XVAL)%in%rownames(coef(lassomodMCV,s=optimal.lambda))]
	CCM <- c(CC,survConcordance(Surv(eval(parse(text=paste0("AMV",k,"$surv"))),eval(parse(text=paste0("AMV",k,"$out"))))~predict(lassomodMCV,newx=XVAL,s=optimal.lambda))$concordance)
	
	CCMT <- c(CCT,survConcordance(Surv(eval(parse(text=paste0("AMV",k,"$surv"))),eval(parse(text=paste0("AMV",k,"$out"))))~predict(lassomodMCV,newx=XVAL,s="lambda.min"))$concordance)
	
	NNM <- c(NN,length(selectedBeta))
	NNMT <- c(NNT,length(selectedBetaT))
}


########################################
## VARIABLE SELECTION FOR FINAL SCORE ##
########################################


# These should be excluded because from verbal interview
to_ex <- c("f.134.0.0*age","f.135.0.0*age","f.137.0.0*age","f.189.0.0")
to_selM <- to_sel[!to_sel%in%to_ex]
formufinM2 <- paste0("~",paste0(to_selM, collapse="+"))

set.seed(123)
MODSTEPW <- cph(as.formula(paste0("Surv(surv,out==1)",formufinM2)),data=AMT1,x=T,y=T)
STEPWM <- fastbw(MODSTEPW,type="residual")

ph <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)
ph[as.character(ph[,1])%in%unlist(lapply(strsplit(STEPWM$names.kept,"\\."),"[",2)),c(6,7)]


# Add extra variables for conformity with uk biobank # 
newaddM <- c(STEPWM$names.kept,FFM[NAMESM %in% c("f.1239.0.0","f.2443.0.0","f.6150.0.l1","f.6150.0.l2","f.6150.0.l3","f.6150.0.l4","f.6150.0.l5","f.6145.0.l1","f.6145.0.l3","f.6145.0.l4","f.6145.0.l5","f.6145.0.l6","f.6145.0.l7","f.6146.0.l2","f.6146.0.l3","f.6146.0.l4","f.709.0.0","f.6141.0.l2","f.6141.0.l3","f.6141.0.l4","f.6141.0.l5","f.6141.0.l6","f.6141.0.l7","f.6141.0.l8")])


# Replace the two questions about Diabetes and heart attack yes/no because self-reported age is tricky #
newaddM <- newaddM[!newaddM%in% c("f.2976.0.0","f.3894.0.0")]

ph[as.character(ph[,1])%in%unlist(lapply(strsplit(newaddM,"\\."),"[",2)),c(1,6,7)]

MAINMODM <- coxph.mids(as.formula(paste0("Surv(surv,out==1)~",paste0(newaddM,collapse="+"))), 5, N_MISSM,"AMT")
write.csv(summary(pool(MAINMODM)), file="test.csv")

save(MAINMODM, file="/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")



#################################################
## GEOGRAPHICAL DISCRIMINATION AND CALIBRATION ## 
#################################################


# DISCRIMINATION #
formufinM3 <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]])

class.temp1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM3)), AMT1)			
class.temp2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM3)), AMT2)
class.temp3 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM3)), AMT3)
class.temp4 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM3)), AMT4)
class.temp5 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM3)), AMT5)

c1 <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1,AMV1))$concordance
c2 <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2,AMV2))$concordance
c3 <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3,AMV3))$concordance
c4 <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4,AMV4))$concordance
c5 <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5,AMV5))$concordance

CINDEXM <- mean(c(c1,c2,c3,c4,c5))


c1S <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1,AMV1))$std.err
c2S <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2,AMV2))$std.err
c3S <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3,AMV3))$std.err
c4S <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4,AMV4))$std.err
c5S <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5,AMV5))$std.err

CINDEXMSE <- mean(c(c1S,c2S,c3S,c4S,c5S))

mean(CINDEXM)+1.96*CINDEXMSE
mean(CINDEXM)-1.96*CINDEXMSE


# CALIBRATION #
Class_riskM1 <- as.numeric(1-summary(survfit(class.temp1,AMV1,se.fit=F), times=5)$surv)
Class_riskM2 <- as.numeric(1-summary(survfit(class.temp2,AMV2,se.fit=F), times=5)$surv)
Class_riskM3 <- as.numeric(1-summary(survfit(class.temp3,AMV3,se.fit=F), times=5)$surv)
Class_riskM4 <- as.numeric(1-summary(survfit(class.temp4,AMV4,se.fit=F), times=5)$surv)
Class_riskM5 <- as.numeric(1-summary(survfit(class.temp5,AMV5,se.fit=F), times=5)$surv)

				
loe1C <- lowess(Class_riskM1, AMV1$out5, iter = 0)
loe2C <- lowess(Class_riskM2, AMV2$out5, iter = 0)
loe3C <- lowess(Class_riskM3, AMV3$out5, iter = 0)
loe4C <- lowess(Class_riskM4, AMV4$out5, iter = 0)
loe5C <- lowess(Class_riskM5, AMV5$out5, iter = 0)

hl1C <- hoslem.test(AMV1$out5, Class_riskM1, g=10)
hlexp1C <- cbind(hl1C$observed[,2]/hl1C$observed[,1],hl1C$expected[,2]/hl1C$expected[,1])

hl2C <- hoslem.test(AMV2$out5, Class_riskM2, g=10)
hlexp2C <- cbind(hl2C$observed[,2]/hl2C$observed[,1],hl2C$expected[,2]/hl2C$expected[,1])

hl3C <- hoslem.test(AMV3$out5, Class_riskM3, g=10)
hlexp3C <- cbind(hl3C$observed[,2]/hl3C$observed[,1],hl3C$expected[,2]/hl3C$expected[,1])

hl4C <- hoslem.test(AMV4$out5, Class_riskM4, g=10)
hlexp4C <- cbind(hl4C$observed[,2]/hl4C$observed[,1],hl4C$expected[,2]/hl4C$expected[,1])

hl5C <- hoslem.test(AMV5$out5, Class_riskM5, g=10)
hlexp5C <- cbind(hl5C$observed[,2]/hl5C$observed[,1],hl5C$expected[,2]/hl5C$expected[,1])


loeEM <- cbind(rowMeans(cbind(loe1C$x,loe2C$x,loe3C$x,loe4C$x,loe5C$x)),rowMeans(cbind(loe1C$y,loe2C$y,loe3C$y,loe4C$y,loe5C$y)))
hlexpEM <- cbind(rowMeans(cbind(hlexp1C[,1],hlexp2C[,1],hlexp3C[,1],hlexp4C[,1],hlexp5C[,1])),rowMeans(cbind(hlexp1C[,2],hlexp2C[,2],hlexp3C[,2],hlexp4C[,2],hlexp5C[,2])))


HLTESTEM <- chisqMI(c(hl1C$statistic,hl2C$statistic,hl3C$statistic,hl4C$statistic,hl5C$statistic),df=8)

save(loeEM,file="/proj/b2011036/uk.biobank/prediction_results/loeEM.Rdata")
save(hlexpEM,file="/proj/b2011036/uk.biobank/prediction_results/hlexpEM.Rdata")



##############
## ONLY AGE ## 
##############

class.temp1AGE <- coxph(Surv(surv,out==1)~age, AMT1)			
c1AGEM <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1AGE,AMV1))$concordance
c1AGEMSE <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1AGE,AMV1))$std.err
c1AGEM+1.96*c1AGEMSE
c1AGEM-1.96*c1AGEMSE

#####################################
## IMPROVEMENT ABOVE CHARSON INDEX ## 
#####################################

formufinM4 <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]],"+charlsonSR")
formufinM5 <- paste0("~","age+charlsonSR")


## Age + Charlson score ##
class.temp1C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT1)			
class.temp2C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT2)
class.temp3C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT3)
class.temp4C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT4)
class.temp5C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT5)

c1C1 <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1C1,AMV1))$concordance
c2C1 <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2C1,AMV2))$concordance
c3C1 <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3C1,AMV3))$concordance
c4C1 <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4C1,AMV4))$concordance
c5C1 <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5C1,AMV5))$concordance

CINDEXMC1 <- mean(c(c1C1,c2C1,c3C1,c4C1,c5C1))

c1C1S <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1C1,AMV1))$std.err
c2C1S <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2C1,AMV2))$std.err
c3C1S <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3C1,AMV3))$std.err
c4C1S <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4C1,AMV4))$std.err
c5C1S <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5C1,AMV5))$std.err

CINDEXMC1SE <- mean(c(c1C1S,c2C1S,c3C1S,c4C1S,c5C1S))

mean(CINDEXMC1)+1.96*CINDEXMC1SE
mean(CINDEXMC1)-1.96*CINDEXMC1SE

## Predictioon score + Charlson score ##
class.temp1C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM4)), AMT1)			
class.temp2C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM4)), AMT2)
class.temp3C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM4)), AMT3)
class.temp4C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM4)), AMT4)
class.temp5C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM4)), AMT5)

c1C2 <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1C2,AMV1))$concordance
c2C2 <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2C2,AMV2))$concordance
c3C2 <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3C2,AMV3))$concordance
c4C2 <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4C2,AMV4))$concordance
c5C2 <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5C2,AMV5))$concordance

CINDEXMC2 <- mean(c(c1C2,c2C2,c3C2,c4C2,c5C2))




#############
### FEMALE ##
#############


NAMESF <- colnames(AF1)[!colnames(AF1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

resF <- read.xls("/proj/b2011036/uk.biobank/univariateF.xlsx", header=T, stringsAsFactor=F)

##EXCLUDE ALL CATEGORIES SHOULD NOT BE INCLUDED IN THE SCORE ##
resF2 <- resF[!resF$Measurement.class%in%c("Physical measures","Blood assays","Cognitive function"),]


## EXTRACT TOP 20 MEASUREMENTS FOR EACH CAUSE OF DEATH
FFF20 <- NULL
for (dis in c("CA","CVD","RE","DG","EC","OD"))
{
	K <- NULL
	for (f in unique(resF2$Original.code))
	{
		
		t <- resF2[resF2$Original.code == f,]
		if(unique(t$P.value.for.Schoenfeld.residuals)!= "&lt;0.00001")
		{
			t2 <- t[t$Age.category=="All",]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),unique(t$Original.code),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}	
					else {k <- c(unique(t$Original.code),unique(t$Original.code),NA)}
							
		}	
		else if(unique(t$P.value.for.Schoenfeld.residuals)== "&lt;0.00001")
		{
			t2 <- t[!t$Age.category%in%c("All","Reference"),]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}		
					else {k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),NA)}		
		}
		
		K <- rbind(K,k)
	}	
	
	FFF20 <- cbind(FFF20,K[order(-as.numeric(K[,3])),2][1:20])
}



to_sel <- unique(as.vector(FFF20))
formufinF <- paste0("~",paste0(to_sel, collapse="+"))



## FORMULA INCLUDING ALL THE FACTORS (THIS IS USED LATER) ##
FFF <- NULL
for (f in NAMESF)
{		
	if (!any(resF$C.index.all[resF$Original.code == f] == "Singular fit") & !any(resF$C.index.all[resF$Original.code == f] == "N <= 5"))
	{
			k <- ifelse(unique(resF$P.value.for.Schoenfeld.residuals[resF$Original.code == f]) == "&lt;0.00001" , paste0(unique(resF$Original.code[resF$Original.code == f]),"*age"), unique(resF$Original.code[resF$Original.code == f]))
	}
	FFF <- c(FFF,k)
}


###############################################################
## LASSO REGRESSION - RUN PARALLELIZED TO SPEED CALCULATIONS ##
###############################################################

# Formula for lasso, including all the measurements
formufinFL <- paste0("~",paste0(FFF, collapse="+"))

registerDoMC(cores=10)
CCF <- NULL
NNF <- NULL
CCFT <- NULL
NNFT <- NULL	
for(k in 1:5)
{
	print(k)
	X <- model.matrix(as.formula(formufinFL), eval(parse(text=paste0("AFT",k))))[,-1]

	set.seed(123)
	# Run penalized cox
	p.fac = rep(1, ncol(X))
	p.fac[which(colnames(X)=="age")] = 0
	
	lassomodFCV <- cv.glmnet(x=X, y=Surv(eval(parse(text=paste0("AFT",k,"$surv"))),eval(parse(text=paste0("AFT",k,"$out")))),family="cox", standardize=F,alpha=1, type.measure="deviance", nfolds=10, penalty.factor = p.fac, nlambda=200, parallel=TRUE)

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
	XVAL <- model.matrix(as.formula(formufinFL), data=eval(parse(text=paste0("AFV",k))))
	XVAL <- XVAL[,colnames(XVAL)%in%rownames(coef(lassomodFCV,s=optimal.lambda))]
	CCF <- c(CCF,survConcordance(Surv(eval(parse(text=paste0("AFV",k,"$surv"))),eval(parse(text=paste0("AFV",k,"$out"))))~predict(lassomodFCV,newx=XVAL,s=optimal.lambda))$concordance)
	
	CCFT <- c(CCFT,survConcordance(Surv(eval(parse(text=paste0("AFV",k,"$surv"))),eval(parse(text=paste0("AFV",k,"$out"))))~predict(lassomodFCV,newx=XVAL,s="lambda.min"))$concordance)
	
	NNF <- c(NN,length(selectedBeta))
	NNFT <- c(NNT,length(selectedBetaT))
}


########################################
## VARIABLE SELECTION FOR FINAL SCORE ##
########################################


# These should be excluded because not measurables
to_ex <- c("f.46.0.0","f.47.0.0","f.102.0.0","f.134.0.0*age","f.135.0.0","f.136.0.0","f.137.0.0")

to_selF <- to_sel[!to_sel%in%to_ex]
formufinF2 <- paste0("~",paste0(to_selF, collapse="+"))

MODSTEPW <- cph(as.formula(paste0("Surv(surv,out==1)",formufinF2)),data=AFT1,x=T,y=T)

set.seed(123)
STEPWF <- fastbw(MODSTEPW,type="residual") 

ph <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)
ph[as.character(ph[,1])%in%unlist(lapply(strsplit(STEPWF$names.kept,"\\."),"[",2)),c(6,7)]


# These variables are added for consistency with UK biobank questions and consistency with men #
newaddF <- c(STEPWF$names.kept,FFF[NAMESF %in% c("f.1239.0.0","f.6145.0.l1","f.6145.0.l3","f.6145.0.l4","f.6145.0.l5","f.6145.0.l6","f.6145.0.l7","f.6146.0.l1","f.6146.0.l2","f.6146.0.l3","f.6146.0.l4")])


# Exclude this question because difficult requires too many extra questions #
newaddF <- newaddF[!newaddF%in% c("f.3456.0.0")]

ph[as.character(ph[,1])%in%unlist(lapply(strsplit(newaddF,"\\."),"[",2)),c(1,6,7)]

MAINMODF <- coxph.mids(as.formula(paste0("Surv(surv,out==1)~",paste0(newaddF,collapse="+"))), 5, N_MISSF,"AFT")
write.csv(summary(pool(MAINMODF)), file="test.csv")

save(MAINMODF, file="/proj/b2011036/uk.biobank/prediction_results/MAINMODF.Rdata")

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODF.Rdata")


#################################################
## GEOGRAPHICAL DISCRIMINATION AND CALIBRATION ## 
#################################################

# DISCRIMINATION #
formufinF3 <- paste0("~",strsplit(as.character(MAINMODF$analyses[[1]]$formula),"~")[[3]])

class.temp1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT1)			
class.temp2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT2)
class.temp3 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT3)
class.temp4 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT4)
class.temp5 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT5)

c1 <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1,AFV1))$concordance
c2 <- survConcordance(Surv(AFV2$surv,AFV2$out==1)~predict(class.temp2,AFV2))$concordance
c3 <- survConcordance(Surv(AFV3$surv,AFV3$out==1)~predict(class.temp3,AFV3))$concordance
c4 <- survConcordance(Surv(AFV4$surv,AFV4$out==1)~predict(class.temp4,AFV4))$concordance
c5 <- survConcordance(Surv(AFV5$surv,AFV5$out==1)~predict(class.temp5,AFV5))$concordance

CINDEXF <- mean(c(c1,c2,c3,c4,c5))


c1S <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1,AFV1))$std.err
c2S <- survConcordance(Surv(AFV2$surv,AFV2$out==1)~predict(class.temp2,AFV2))$std.err
c3S <- survConcordance(Surv(AFV3$surv,AFV3$out==1)~predict(class.temp3,AFV3))$std.err
c4S <- survConcordance(Surv(AFV4$surv,AFV4$out==1)~predict(class.temp4,AFV4))$std.err
c5S <- survConcordance(Surv(AFV5$surv,AFV5$out==1)~predict(class.temp5,AFV5))$std.err

CINDEXFSE <- mean(c(c1S,c2S,c3S,c4S,c5S))

mean(CINDEXF)+1.96*CINDEXFSE
mean(CINDEXF)-1.96*CINDEXFSE


# CALIBRATION #
Class_riskF1 <- as.numeric(1-summary(survfit(class.temp1,AFV1,se.fit=F), times=5)$surv)
Class_riskF2 <- as.numeric(1-summary(survfit(class.temp2,AFV2,se.fit=F), times=5)$surv)
Class_riskF3 <- as.numeric(1-summary(survfit(class.temp3,AFV3,se.fit=F), times=5)$surv)
Class_riskF4 <- as.numeric(1-summary(survfit(class.temp4,AFV4,se.fit=F), times=5)$surv)
Class_riskF5 <- as.numeric(1-summary(survfit(class.temp5,AFV5,se.fit=F), times=5)$surv)

				
loe1C <- lowess(Class_riskF1, AFV1$out5, iter = 0)
loe2C <- lowess(Class_riskF2, AFV2$out5, iter = 0)
loe3C <- lowess(Class_riskF3, AFV3$out5, iter = 0)
loe4C <- lowess(Class_riskF4, AFV4$out5, iter = 0)
loe5C <- lowess(Class_riskF5, AFV5$out5, iter = 0)

hl1C <- hoslem.test(AFV1$out5, Class_riskF1, g=10)
hlexp1C <- cbind(hl1C$observed[,2]/hl1C$observed[,1],hl1C$expected[,2]/hl1C$expected[,1])

hl2C <- hoslem.test(AFV2$out5, Class_riskF2, g=10)
hlexp2C <- cbind(hl2C$observed[,2]/hl2C$observed[,1],hl2C$expected[,2]/hl2C$expected[,1])

hl3C <- hoslem.test(AFV3$out5, Class_riskF3, g=10)
hlexp3C <- cbind(hl3C$observed[,2]/hl3C$observed[,1],hl3C$expected[,2]/hl3C$expected[,1])

hl4C <- hoslem.test(AFV4$out5, Class_riskF4, g=10)
hlexp4C <- cbind(hl4C$observed[,2]/hl4C$observed[,1],hl4C$expected[,2]/hl4C$expected[,1])

hl5C <- hoslem.test(AFV5$out5, Class_riskF5, g=10)
hlexp5C <- cbind(hl5C$observed[,2]/hl5C$observed[,1],hl5C$expected[,2]/hl5C$expected[,1])


loeEF <- cbind(rowMeans(cbind(loe1C$x,loe2C$x,loe3C$x,loe4C$x,loe5C$x)),rowMeans(cbind(loe1C$y,loe2C$y,loe3C$y,loe4C$y,loe5C$y)))
hlexpEF <- cbind(rowMeans(cbind(hlexp1C[,1],hlexp2C[,1],hlexp3C[,1],hlexp4C[,1],hlexp5C[,1])),rowMeans(cbind(hlexp1C[,2],hlexp2C[,2],hlexp3C[,2],hlexp4C[,2],hlexp5C[,2])))


HLTESTEF <- chisqMI(c(hl1C$statistic,hl2C$statistic,hl3C$statistic,hl4C$statistic,hl5C$statistic),df=8)

save(loeEF,file="/proj/b2011036/uk.biobank/prediction_results/loeEF.Rdata")
save(hlexpEF,file="/proj/b2011036/uk.biobank/prediction_results/hlexpEF.Rdata")


##############
## ONLY AGE ## 
##############

class.temp1AGE <- coxph(Surv(surv,out==1)~age, AFT1)			
c1AGEF <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1AGE,AFV1))$concordance
c1AGEFSE <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1AGE,AFV1))$std.err
c1AGEF+1.96*c1AGEFSE
c1AGEF-1.96*c1AGEFSE


#####################################
## IMPROVEMENT ABOVE CHARSON INDEX ## 
#####################################

formufinF4 <- paste0("~",strsplit(as.character(MAINMODF$analyses[[1]]$formula),"~")[[3]],"+charlsonSR")
formufinF5 <- paste0("~","age+charlsonSR")


## Age + Charlson score ##
class.temp1C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT1)			
class.temp2C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT2)
class.temp3C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT3)
class.temp4C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT4)
class.temp5C1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT5)

c1C1 <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1C1,AFV1))$concordance
c2C1 <- survConcordance(Surv(AFV2$surv,AFV2$out==1)~predict(class.temp2C1,AFV2))$concordance
c3C1 <- survConcordance(Surv(AFV3$surv,AFV3$out==1)~predict(class.temp3C1,AFV3))$concordance
c4C1 <- survConcordance(Surv(AFV4$surv,AFV4$out==1)~predict(class.temp4C1,AFV4))$concordance
c5C1 <- survConcordance(Surv(AFV5$surv,AFV5$out==1)~predict(class.temp5C1,AFV5))$concordance

CINDEXFC1 <- mean(c(c1C1,c2C1,c3C1,c4C1,c5C1))

c1C1S <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1C1,AFV1))$std.err
c2C1S <- survConcordance(Surv(AFV2$surv,AFV2$out==1)~predict(class.temp2C1,AFV2))$std.err
c3C1S <- survConcordance(Surv(AFV3$surv,AFV3$out==1)~predict(class.temp3C1,AFV3))$std.err
c4C1S <- survConcordance(Surv(AFV4$surv,AFV4$out==1)~predict(class.temp4C1,AFV4))$std.err
c5C1S <- survConcordance(Surv(AFV5$surv,AFV5$out==1)~predict(class.temp5C1,AFV5))$std.err

CINDEXFC1SE <- mean(c(c1C1S,c2C1S,c3C1S,c4C1S,c5C1S))

mean(CINDEXFC1)+1.96*CINDEXFC1SE
mean(CINDEXFC1)-1.96*CINDEXFC1SE


## Predictioon score + Charlson score ##
class.temp1C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF4)), AFT1)			
class.temp2C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF4)), AFT2)
class.temp3C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF4)), AFT3)
class.temp4C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF4)), AFT4)
class.temp5C2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF4)), AFT5)

c1C2 <- survConcordance(Surv(AFV1$surv,AFV1$out==1)~predict(class.temp1C2,AFV1))$concordance
c2C2 <- survConcordance(Surv(AFV2$surv,AFV2$out==1)~predict(class.temp2C2,AFV2))$concordance
c3C2 <- survConcordance(Surv(AFV3$surv,AFV3$out==1)~predict(class.temp3C2,AFV3))$concordance
c4C2 <- survConcordance(Surv(AFV4$surv,AFV4$out==1)~predict(class.temp4C2,AFV4))$concordance
c5C2 <- survConcordance(Surv(AFV5$surv,AFV5$out==1)~predict(class.temp5C2,AFV5))$concordance

CINDEXFC2 <- mean(c(c1C2,c2C2,c3C2,c4C2,c5C2))



save(loeEF,file="/proj/b2011036/uk.biobank/prediction_results/loeEF.Rdata")
save(hlexpEF,file="/proj/b2011036/uk.biobank/prediction_results/hlexpEF.Rdata")




#######################################################################
#######################################################################
####### 6. RE-WEIGHT BY MORTALITY TABLES AND CENSUS INFORMATION #######
#######################################################################
#######################################################################



######################################################
### WEIGHTS TO REWEIGHT THE COVARIATE DISTRIBUTION  ##
######################################################

GHM <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/general_health_M.csv", stringsAsFactor=F) 
GHF <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/general_health_F.csv", stringsAsFactor=F) 
VCM <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/van_car_M.csv", stringsAsFactor=F) 
VCF <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/van_car_F.csv", stringsAsFactor=F) 

### SELF REPORTED HEALTH ###

GHMp <- NULL
GHFp <- NULL
for (i in 1:nrow(GHM))
{
	GHMp <- rbind(GHMp, c(GHM[i,2]/sum(GHM[i,2:5]),GHM[i,3]/sum(GHM[i,2:5]),GHM[i,4]/sum(GHM[i,2:5]),GHM[i,5]/sum(GHM[i,2:5])))
	GHFp <- rbind(GHFp, c(GHF[i,2]/sum(GHF[i,2:5]),GHF[i,3]/sum(GHF[i,2:5]),GHF[i,4]/sum(GHF[i,2:5]),GHM[i,5]/sum(GHF[i,2:5])))
}


# Get the weights ##
ttM <- table(AMT1$f.2178.0.0,AMT1$age)
ttF <- table(AFT1$f.2178.0.0,AFT1$age)

GHUKMp <- t(ttM)/rowSums(t(ttM))
GHUKFp <- t(ttF)/rowSums(t(ttF))


# Calculate the new weighted average to use as model mean in the assessment of risk #
GHWM <- GHMp/GHUKMp
GHMnewmean <- colSums(GHWM*t(ttM))/sum(colSums(GHWM*t(ttM)))

GHWF <- GHFp/GHUKFp
GHFnewmean <- colSums(GHWF*t(ttF))/sum(colSums(GHWF*t(ttF)))

attr(GHMnewmean,"name") <- "f.2178.0.0"
save(GHMnewmean,file="/proj/b2011036/uk.biobank/out_weights/GHMnewmean.Rdata")
attr(GHFnewmean,"name") <- "f.2178.0.0"
save(GHFnewmean,file="/proj/b2011036/uk.biobank/out_weights/GHFnewmean.Rdata")



### NUMBER OF CARS OR VAN ###

VCMp <- NULL
VCFp <- NULL
for (i in 1:nrow(VCM))
{
	VCMp <- rbind(VCMp, c(VCM[i,2]/sum(VCM[i,2:6]),VCM[i,3]/sum(VCM[i,2:6]),VCM[i,4]/sum(VCM[i,2:6]),VCM[i,5]/sum(VCM[i,2:6]),VCM[i,6]/sum(VCM[i,2:6])))
	VCFp <- rbind(VCFp, c(VCF[i,2]/sum(VCF[i,2:6]),VCF[i,3]/sum(VCF[i,2:6]),VCF[i,4]/sum(VCF[i,2:6]),VCF[i,5]/sum(VCF[i,2:6]),VCF[i,6]/sum(VCF[i,2:6])))
}



# Distribution in the Uk biobank
ttM <- table(AMT1$f.728.0.0,AMT1$age)
ttF <- table(AFT1$f.728.0.0,AFT1$age)


N_vcM <- rbind(ttM[1,]+ttM[6,],ttM[2,],ttM[3,],ttM[4,],ttM[5,])
rownames(N_vcM) <- c("0","1","2","3","4+")
VCUKMp <- t(N_vcM)/rowSums(t(N_vcM))

N_vcF <- rbind(ttF[1,]+ttF[6,],ttF[2,],ttF[3,],ttF[4,],ttF[5,])
rownames(N_vcF) <- c("0","1","2","3","4+")
VCUKFp <- t(N_vcF)/rowSums(t(N_vcF))


# Calculate the new weighted average to use as model mean in the assessment of risk #
VCWM <- VCMp/VCUKMp
VCWM <- cbind(VCWM[,1],VCWM[,2],VCWM[,3],VCWM[,4],VCWM[,5],VCWM[,1])
VCMnewmean <- colSums(VCWM*t(ttM))/sum(colSums(VCWM*t(ttM)))

VCWF <- VCFp/VCUKFp
VCWF <- cbind(VCWF[,1],VCWF[,2],VCWF[,3],VCWF[,4],VCWF[,5],VCWF[,1])
VCFnewmean <- colSums(VCWF*t(ttF))/sum(colSums(VCWF*t(ttF)))

# Remove no-home, these are not considered
VCMnewmean <- c(VCMnewmean[1],VCMnewmean[2],VCMnewmean[3],VCMnewmean[4],VCMnewmean[5])
attr(VCMnewmean,"name") <- "f.728.0.0"

VCFnewmean <- c(VCFnewmean[1],VCFnewmean[2],VCFnewmean[3],VCFnewmean[4],VCFnewmean[5])
attr(VCFnewmean,"name") <- "f.728.0.0"


attr(VCMnewmean,"name") <- "f.728.0.0"
save(VCMnewmean,file="/proj/b2011036/uk.biobank/out_weights/VCMnewmean.Rdata")
attr(VCFnewmean,"name") <- "f.728.0.0"
save(VCFnewmean,file="/proj/b2011036/uk.biobank/out_weights/VCFnewmean.Rdata")


#######################################
#### WEIGHTS FROM MORTALITY TABLES ####
#######################################


#################################
### OVERALL MORTALITY WEIGHTS ###
#################################

muktM <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/male_2009_11.csv")
muktF <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/female_2009_11.csv")
maxsurv <- 5
# Male #
MM <- NULL
FF <- NULL
SF5age <- NULL
SM5age <- NULL
# KM, time-in-study
kmMage <- survfit(Surv(age,age_exit, out)~1,data=AM1,se.fit=FALSE)
kmFage <- survfit(Surv(age,age_exit, out)~1,data=AF1,se.fit=FALSE)

for (i in 40:70)
{
	# Uk lifetables
	M <- muktM[muktM[,1]%in%(i+maxsurv),4]/muktM[muktM[,1]%in%i,4]
	F <- muktF[muktF[,1]%in%(i+maxsurv),4]/muktF[muktF[,1]%in%i,4]
	# KM age as time-scale
	SM5age <- c(SM5age, kmMage$surv[which(abs(kmMage$time-(i+5))==min(abs(kmMage$time-(i+5))))]/kmMage$surv[which(abs(kmMage$time-i)==min(abs(kmMage$time-i)))])
	SF5age <- c(SF5age, kmFage$surv[which(abs(kmFage$time-(i+5))==min(abs(kmFage$time-(i+5))))]/kmFage$surv[which(abs(kmFage$time-i)==min(abs(kmFage$time-i)))])
	MM <- c(MM,M)
	FF <- c(FF,F)
}

wwMage <- -log(MM)/-log(SM5age)
names(wwMage) <- 40:70
wwFage <- -log(FF)/-log(SF5age)
names(wwFage) <- 40:70

save(wwMage, file="/proj/b2011036/uk.biobank/out_weights/wwMage.Rdata")
save(wwFage, file="/proj/b2011036/uk.biobank/out_weights/wwFage.Rdata")


## EXPORT QUANTITY FOR SHINY ##
#  5-YEARS SURVIVAL USING UK TABLE (FROM 0 to 95 YEARS) #
FF095 <- NULL
MM095 <- NULL
for (i in 15:95)
{
M <- muktM[muktM[,1]%in%(i+maxsurv),4]/muktM[muktM[,1]%in%i,4]
F <- muktF[muktF[,1]%in%(i+maxsurv),4]/muktF[muktF[,1]%in%i,4]
MM095 <- c(MM095,M)
FF095 <- c(FF095,F)
}
names(MM095) <- 15:95
names(FF095) <- 15:95

save(MM095, file="/proj/b2011036/uk.biobank/out_weights/MM095.Rdata")
save(FF095, file="/proj/b2011036/uk.biobank/out_weights/FF095.Rdata")





#######################################################################
#######################################################################
####### 7. CALCULATE THE RE-CALIBRATED INDIVIDUAL RISK ################
####### AND SAVE MAIN QUANTITIES FOR PYTHON SCRIPT ####################
#######################################################################
#######################################################################


##############
#### MALE ####
##############

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

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")

load("/proj/b2011036/uk.biobank/out_weights/VCMnewmean.Rdata")
load("/proj/b2011036/uk.biobank/out_weights/GHMnewmean.Rdata")
load("/proj/b2011036/uk.biobank/out_weights/wwMage.Rdata")



## RECODE THE DATASET TO:
# 1. EXCLUDE THOSE NOT WALKING OR REPLYING NONE OF THE ABOVE TO QUESTION 924

AMT1EX <- AMT1[!AMT1$f.924.0.0%in%c("Not walking") & AMT1$f.709.0.0!="No home" & AMT1$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMT1EX$f.924.0.0 <- relevel(AMT1EX$f.924.0.0, ref="Slow pace")
AMT1EX$f.924.0.0 <- factor(AMT1EX$f.924.0.0)
AMT1EX$f.709.0.0 <- relevel(AMT1EX$f.709.0.0, ref="[1,2]")
AMT1EX$f.709.0.0 <- factor(AMT1EX$f.709.0.0)
AMT1EX$f.728.0.0 <- relevel(AMT1EX$f.728.0.0, ref="None")
AMT1EX$f.728.0.0 <- factor(AMT1EX$f.728.0.0)


#Do the same for validation
AMV1EX <- AMV1[!AMV1$f.924.0.0%in%c("Not walking") & AMV1$f.709.0.0!="No home" & AMV1$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMV1EX$f.924.0.0 <- relevel(AMV1EX$f.924.0.0, ref="Slow pace")
AMV1EX$f.924.0.0 <- factor(AMV1EX$f.924.0.0)
AMV1EX$f.709.0.0 <- relevel(AMV1EX$f.709.0.0, ref="[1,2]")
AMV1EX$f.709.0.0 <- factor(AMV1EX$f.709.0.0)
AMV1EX$f.728.0.0 <- relevel(AMV1EX$f.728.0.0, ref="None")
AMV1EX$f.728.0.0 <- factor(AMV1EX$f.728.0.0)


AMT2EX <- AMT2[!AMT2$f.924.0.0%in%c("Not walking") & AMT2$f.709.0.0!="No home" & AMT2$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMT2EX$f.924.0.0 <- relevel(AMT2EX$f.924.0.0, ref="Slow pace")
AMT2EX$f.924.0.0 <- factor(AMT2EX$f.924.0.0)
AMT2EX$f.709.0.0 <- relevel(AMT2EX$f.709.0.0, ref="[1,2]")
AMT2EX$f.709.0.0 <- factor(AMT2EX$f.709.0.0)
AMT2EX$f.728.0.0 <- relevel(AMT2EX$f.728.0.0, ref="None")
AMT2EX$f.728.0.0 <- factor(AMT2EX$f.728.0.0)


#Do the same for validation
AMV2EX <- AMV2[!AMV2$f.924.0.0%in%c("Not walking") & AMV2$f.709.0.0!="No home" & AMV2$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMV2EX$f.924.0.0 <- relevel(AMV2EX$f.924.0.0, ref="Slow pace")
AMV2EX$f.924.0.0 <- factor(AMV2EX$f.924.0.0)
AMV2EX$f.709.0.0 <- relevel(AMV2EX$f.709.0.0, ref="[1,2]")
AMV2EX$f.709.0.0 <- factor(AMV2EX$f.709.0.0)
AMV2EX$f.728.0.0 <- relevel(AMV2EX$f.728.0.0, ref="None")
AMV2EX$f.728.0.0 <- factor(AMV2EX$f.728.0.0)


AMT3EX <- AMT3[!AMT3$f.924.0.0%in%c("Not walking") & AMT3$f.709.0.0!="No home" & AMT3$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMT3EX$f.924.0.0 <- relevel(AMT3EX$f.924.0.0, ref="Slow pace")
AMT3EX$f.924.0.0 <- factor(AMT3EX$f.924.0.0)
AMT3EX$f.709.0.0 <- relevel(AMT3EX$f.709.0.0, ref="[1,2]")
AMT3EX$f.709.0.0 <- factor(AMT3EX$f.709.0.0)
AMT3EX$f.728.0.0 <- relevel(AMT3EX$f.728.0.0, ref="None")
AMT3EX$f.728.0.0 <- factor(AMT3EX$f.728.0.0)


#Do the same for validation
AMV3EX <- AMV3[!AMV3$f.924.0.0%in%c("Not walking") & AMV3$f.709.0.0!="No home" & AMV3$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMV3EX$f.924.0.0 <- relevel(AMV3EX$f.924.0.0, ref="Slow pace")
AMV3EX$f.924.0.0 <- factor(AMV3EX$f.924.0.0)
AMV3EX$f.709.0.0 <- relevel(AMV3EX$f.709.0.0, ref="[1,2]")
AMV3EX$f.709.0.0 <- factor(AMV3EX$f.709.0.0)
AMV3EX$f.728.0.0 <- relevel(AMV3EX$f.728.0.0, ref="None")
AMV3EX$f.728.0.0 <- factor(AMV3EX$f.728.0.0)


AMT4EX <- AMT4[!AMT4$f.924.0.0%in%c("Not walking") & AMT4$f.709.0.0!="No home" & AMT4$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMT4EX$f.924.0.0 <- relevel(AMT4EX$f.924.0.0, ref="Slow pace")
AMT4EX$f.924.0.0 <- factor(AMT4EX$f.924.0.0)
AMT4EX$f.709.0.0 <- relevel(AMT4EX$f.709.0.0, ref="[1,2]")
AMT4EX$f.709.0.0 <- factor(AMT4EX$f.709.0.0)
AMT4EX$f.728.0.0 <- relevel(AMT4EX$f.728.0.0, ref="None")
AMT4EX$f.728.0.0 <- factor(AMT4EX$f.728.0.0)


#Do the same for validation
AMV4EX <- AMV4[!AMV4$f.924.0.0%in%c("Not walking") & AMV4$f.709.0.0!="No home" & AMV4$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMV4EX$f.924.0.0 <- relevel(AMV4EX$f.924.0.0, ref="Slow pace")
AMV4EX$f.924.0.0 <- factor(AMV4EX$f.924.0.0)
AMV4EX$f.709.0.0 <- relevel(AMV4EX$f.709.0.0, ref="[1,2]")
AMV4EX$f.709.0.0 <- factor(AMV4EX$f.709.0.0)
AMV4EX$f.728.0.0 <- relevel(AMV4EX$f.728.0.0, ref="None")
AMV4EX$f.728.0.0 <- factor(AMV4EX$f.728.0.0)


AMT5EX <- AMT5[!AMT5$f.924.0.0%in%c("Not walking") & AMT5$f.709.0.0!="No home" & AMT5$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMT5EX$f.924.0.0 <- relevel(AMT5EX$f.924.0.0, ref="Slow pace")
AMT5EX$f.924.0.0 <- factor(AMT5EX$f.924.0.0)
AMT5EX$f.709.0.0 <- relevel(AMT5EX$f.709.0.0, ref="[1,2]")
AMT5EX$f.709.0.0 <- factor(AMT5EX$f.709.0.0)
AMT5EX$f.728.0.0 <- relevel(AMT5EX$f.728.0.0, ref="None")
AMT5EX$f.728.0.0 <- factor(AMT5EX$f.728.0.0)

#Do the same for validation
AMV5EX <- AMV5[!AMV5$f.924.0.0%in%c("Not walking") & AMV5$f.709.0.0!="No home" & AMV5$f.728.0.0!="No home",]

# Make the reference level "slow pace"
AMV5EX$f.924.0.0 <- relevel(AMV5EX$f.924.0.0, ref="Slow pace")
AMV5EX$f.924.0.0 <- factor(AMV5EX$f.924.0.0)
AMV5EX$f.709.0.0 <- relevel(AMV5EX$f.709.0.0, ref="[1,2]")
AMV5EX$f.709.0.0 <- factor(AMV5EX$f.709.0.0)
AMV5EX$f.728.0.0 <- relevel(AMV5EX$f.728.0.0, ref="None")
AMV5EX$f.728.0.0 <- factor(AMV5EX$f.728.0.0)

# Find subject included in all the imputation datasets #
IDcomT <- Reduce(intersect, list(AMT1EX$f.eid,AMT2EX$f.eid,AMT3EX$f.eid,AMT4EX$f.eid,AMT5EX$f.eid))
AMT1EX <- AMT1EX[AMT1EX$f.eid%in%IDcomT,]
AMT2EX <- AMT2EX[AMT2EX$f.eid%in%IDcomT,]
AMT3EX <- AMT3EX[AMT3EX$f.eid%in%IDcomT,]
AMT4EX <- AMT4EX[AMT4EX$f.eid%in%IDcomT,]
AMT5EX <- AMT5EX[AMT5EX$f.eid%in%IDcomT,]

IDcomV <- Reduce(intersect,list(AMV1EX$f.eid,AMV2EX$f.eid,AMV3EX$f.eid,AMV4EX$f.eid,AMV5EX$f.eid))
AMV1EX <- AMV1EX[AMV1EX$f.eid%in%IDcomV,]
AMV2EX <- AMV2EX[AMV2EX$f.eid%in%IDcomV,]
AMV3EX <- AMV3EX[AMV3EX$f.eid%in%IDcomV,]
AMV4EX <- AMV4EX[AMV4EX$f.eid%in%IDcomV,]
AMV5EX <- AMV5EX[AMV5EX$f.eid%in%IDcomV,]


##### EXTRACT COEFFICIENTS AND MEANS ######

# Forumla
formufinM3 <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]])


# Coefficients
coxCC.temp1M <- coxph(as.formula(paste0("Surv(surv,out)",formufinM3)),data=AMT1EX)
coxCC.temp2M <- coxph(as.formula(paste0("Surv(surv,out)",formufinM3)),data=AMT2EX)
coxCC.temp3M <- coxph(as.formula(paste0("Surv(surv,out)",formufinM3)),data=AMT3EX)
coxCC.temp4M <- coxph(as.formula(paste0("Surv(surv,out)",formufinM3)),data=AMT4EX)
coxCC.temp5M <- coxph(as.formula(paste0("Surv(surv,out)",formufinM3)),data=AMT5EX)


ss1 <- summary(coxCC.temp1M)$coefficients[,1]
ss2 <- summary(coxCC.temp2M)$coefficients[,1]
ss3 <- summary(coxCC.temp3M)$coefficients[,1]
ss4 <- summary(coxCC.temp4M)$coefficients[,1]
ss5 <- summary(coxCC.temp5M)$coefficients[,1]

coefM <- rowMeans(cbind(ss1,ss2,ss3,ss4,ss5))

# Means
for (i in 1:length(list(VCMnewmean,GHMnewmean)))
	{coxCC.temp1M$means[grep(attributes(list(VCMnewmean,GHMnewmean)[[i]])$name,names(coxCC.temp1M$means))] <- list(VCMnewmean,GHMnewmean)[[i]][-1]}
mm1 <- coxCC.temp1M$means

for (i in 1:length(list(VCMnewmean,GHMnewmean)))
	{coxCC.temp2M$means[grep(attributes(list(VCMnewmean,GHMnewmean)[[i]])$name,names(coxCC.temp2M$means))] <- list(VCMnewmean,GHMnewmean)[[i]][-1]}
mm2 <- coxCC.temp1M$means

for (i in 1:length(list(VCMnewmean,GHMnewmean)))
	{coxCC.temp3M$means[grep(attributes(list(VCMnewmean,GHMnewmean)[[i]])$name,names(coxCC.temp3M$means))] <- list(VCMnewmean,GHMnewmean)[[i]][-1]}
mm3 <- coxCC.temp1M$means

for (i in 1:length(list(VCMnewmean,GHMnewmean)))
	{coxCC.temp4M$means[grep(attributes(list(VCMnewmean,GHMnewmean)[[i]])$name,names(coxCC.temp4M$means))] <- list(VCMnewmean,GHMnewmean)[[i]][-1]}
mm4 <- coxCC.temp1M$means

for (i in 1:length(list(VCMnewmean,GHMnewmean)))
	{coxCC.temp5M$means[grep(attributes(list(VCMnewmean,GHMnewmean)[[i]])$name,names(coxCC.temp5M$means))] <- list(VCMnewmean,GHMnewmean)[[i]][-1]}
mm5 <- coxCC.temp1M$means

meansM <- rowMeans(cbind(mm1,mm2,mm3,mm4,mm5))


# Variance-covariance matrix - same as pool(model)$ubar
vs1 <- vcov(coxCC.temp1M)
vs2 <- vcov(coxCC.temp2M)
vs3 <- vcov(coxCC.temp3M)
vs4 <- vcov(coxCC.temp4M)
vs5 <- vcov(coxCC.temp5M)

x <- array(c(vs1,vs2,vs3,vs4,vs5), dim=c(72,72,5)) 
vcovM <- apply(x,c(1,2),mean)

# Write data (ordered)
ordetobe <- c("age","f.728.0.0One","f.728.0.0Two","f.728.0.0Three","f.728.0.0Four or more","f.709.0.0(2,4]","f.709.0.0(4,100]","f.6141.0.l12","f.6141.0.l13","f.6141.0.l22","f.6141.0.l23","f.6141.0.l32","f.6141.0.l33","f.6141.0.l42","f.6141.0.l43","f.6141.0.l52","f.6141.0.l53","f.6141.0.l62","f.6141.0.l63","f.6141.0.l72","f.6141.0.l73","f.6141.0.l82","f.6141.0.l83","f.1239.0.02","f.1239.0.03","f.1249.0.02","f.1249.0.03","f.1249.0.04","f.1249.0.05","f.2178.0.02","f.2178.0.03","f.2178.0.04","f.924.0.0Steady average pace","f.924.0.0Brisk pace","f.924.0.0None of the above","f.2443.0.02","f.2453.0.02","f.6150.0.l12","f.6150.0.l22","f.6150.0.l32","f.6150.0.l42","f.6150.0.l52","f.6145.0.l12","f.6145.0.l22","f.6145.0.l32","f.6145.0.l42","f.6145.0.l5Death of a spouse or partner","f.6145.0.l62","f.6145.0.l72","f.6146.0.l12","f.6146.0.l22","f.6146.0.l32","f.6146.0.l42","age:f.709.0.0(2,4]","age:f.709.0.0(4,100]","age:f.6141.0.l22","age:f.6141.0.l23","age:f.6141.0.l42","age:f.6141.0.l43","age:f.6141.0.l62","age:f.6141.0.l63","age:f.6141.0.l82","age:f.6141.0.l83","age:f.924.0.0Steady average pace","age:f.924.0.0Brisk pace","age:f.924.0.0None of the above","age:f.2443.0.02","f.2453.0.02:age","age:f.6150.0.l52","age:f.6146.0.l12","age:f.6146.0.l32","age:f.6146.0.l42")

vcovMO <- vcovM[order(match(colnames(vs1),ordetobe)), order(match(colnames(vs1),ordetobe))]
coefMO <- coefM[order(match(names(coefM),ordetobe))]
meansMO <- meansM[order(match(names(meansM),ordetobe))]


write.csv(vcovMO,file="/proj/b2011036/uk.biobank/out_weights/varcovM.csv")
write.csv(cbind(coefMO,meansMO),file="/proj/b2011036/uk.biobank/out_weights/coefmeanM.csv")


###### BASELINE HAZARD ######
bb1 <- basehaz(coxCC.temp1M)$hazard[which(abs(basehaz(coxCC.temp1M)$time-5)==min(abs(basehaz(coxCC.temp1M)$time-5)))]
bb2 <-basehaz(coxCC.temp2M)$hazard[which(abs(basehaz(coxCC.temp2M)$time-5)==min(abs(basehaz(coxCC.temp2M)$time-5)))]
bb3 <-basehaz(coxCC.temp3M)$hazard[which(abs(basehaz(coxCC.temp3M)$time-5)==min(abs(basehaz(coxCC.temp3M)$time-5)))]
bb4 <-basehaz(coxCC.temp4M)$hazard[which(abs(basehaz(coxCC.temp4M)$time-5)==min(abs(basehaz(coxCC.temp4M)$time-5)))]
bb5 <-basehaz(coxCC.temp5M)$hazard[which(abs(basehaz(coxCC.temp5M)$time-5)==min(abs(basehaz(coxCC.temp5M)$time-5)))]

basehazM <- mean(c(bb1,bb2,bb3,bb4,bb5))

write.csv(basehazM,file="/proj/b2011036/uk.biobank/out_weights/basehazM.csv")


## Calculate the risk in the geographical validation ##

All_riskMT <- mclapply(seq(1:5), 		function(c) { calcweightedRISK(object=eval(parse(text=paste0("coxCC.temp",c,"M"))), newdata = eval(parse(text=paste0("AMV",c,"EX"))),newmeans=list(VCMnewmean,GHMnewmean),weights=wwMage)}, mc.cores=5)

All_riskM <- exp(-exp(rowMeans(cbind(log(-log(All_riskMT[[1]])),
log(-log(All_riskMT[[2]])),
log(-log(All_riskMT[[3]])),
log(-log(All_riskMT[[4]])),
log(-log(All_riskMT[[5]]))))))


save(All_riskM,file="/proj/b2011036/uk.biobank/prediction_results/All_riskM.Rdata")
save(AMV1EX,file="/proj/b2011036/uk.biobank/prediction_results/AMV1EX.Rdata")


## For supplementary table 4 #
ind <- which(AMV1EX$sex=="Male" & AMV1EX$age=="52" & AMV1EX$f.728.0.0=="One" & AMV1EX$f.709.0.0=="[1,2]" &
 AMV1EX$f.6141.0.l1 == "Husband, wife or partner" & AMV1EX$f.6141.0.l2 == "0" & AMV1EX$f.6141.0.l3 == "0" & AMV1EX$f.6141.0.l4 == "0" & AMV1EX$f.6141.0.l5 == "0" & AMV1EX$f.6141.0.l6 == "0" & AMV1EX$f.6141.0.l7 == "0"
& AMV1EX$f.6141.0.l8 == "0" & AMV1EX$f.1239.0.0 == "Yes, on most or all days" & AMV1EX$f.1249.0.0 == "I'm a smoker" & AMV1EX$f.2178.0.0=="Poor" & AMV1EX$f.924.0.0=="Steady average pace" & AMV1EX$f.2443.0.0=="No" & AMV1EX$f.2453.0.0=="No" & AMV1EX$f.6150.0.l1=="None of the above" &  AMV1EX$f.6145.0.l1=="0" & AMV1EX$f.6145.0.l2=="Serious illness, injury or assault to yourself" & AMV1EX$f.6145.0.l3=="0" & AMV1EX$f.6145.0.l4=="Death of a close relative" & AMV1EX$f.6145.0.l5=="0" & AMV1EX$f.6145.0.l6=="0" & AMV1EX$f.6145.0.l7=="Financial difficulties" & AMV1EX$f.6146.0.l1=="None of the above")


All_riskM[ind]



################
#### FEMALE ####
################

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

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODF.Rdata")

load("/proj/b2011036/uk.biobank/out_weights/VCFnewmean.Rdata")
load("/proj/b2011036/uk.biobank/out_weights/GHFnewmean.Rdata")
load("/proj/b2011036/uk.biobank/out_weights/wwFage.Rdata")

## RECODE THE DATASET TO:
# 1. EXCLUDE THOSE NOT WALKING OR REPLYING NONE OF THE ABOVE TO QUESTION 924

AFT1EX <- AFT1[!AFT1$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFT1EX$f.924.0.0 <- relevel(AFT1EX$f.924.0.0, ref="Slow pace")
AFT1EX$f.924.0.0 <- factor(AFT1EX$f.924.0.0)

#Do the same for validation
AFV1EX <- AFV1[!AFV1$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFV1EX$f.924.0.0 <- relevel(AFV1EX$f.924.0.0, ref="Slow pace")
AFV1EX$f.924.0.0 <- factor(AFV1EX$f.924.0.0)


AFT2EX <- AFT2[!AFT2$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFT2EX$f.924.0.0 <- relevel(AFT2EX$f.924.0.0, ref="Slow pace")
AFT2EX$f.924.0.0 <- factor(AFT2EX$f.924.0.0)


#Do the same for validation
AFV2EX <- AFV2[!AFV2$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFV2EX$f.924.0.0 <- relevel(AFV2EX$f.924.0.0, ref="Slow pace")
AFV2EX$f.924.0.0 <- factor(AFV2EX$f.924.0.0)


AFT3EX <- AFT3[!AFT3$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFT3EX$f.924.0.0 <- relevel(AFT3EX$f.924.0.0, ref="Slow pace")
AFT3EX$f.924.0.0 <- factor(AFT3EX$f.924.0.0)

#Do the same for validation
AFV3EX <- AFV3[!AFV3$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFV3EX$f.924.0.0 <- relevel(AFV3EX$f.924.0.0, ref="Slow pace")
AFV3EX$f.924.0.0 <- factor(AFV3EX$f.924.0.0)

AFT4EX <- AFT4[!AFT4$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFT4EX$f.924.0.0 <- relevel(AFT4EX$f.924.0.0, ref="Slow pace")
AFT4EX$f.924.0.0 <- factor(AFT4EX$f.924.0.0)

#Do the same for validation
AFV4EX <- AFV4[!AFV4$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFV4EX$f.924.0.0 <- relevel(AFV4EX$f.924.0.0, ref="Slow pace")
AFV4EX$f.924.0.0 <- factor(AFV4EX$f.924.0.0)


AFT5EX <- AFT5[!AFT5$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFT5EX$f.924.0.0 <- relevel(AFT5EX$f.924.0.0, ref="Slow pace")
AFT5EX$f.924.0.0 <- factor(AFT5EX$f.924.0.0)

#Do the same for validation
AFV5EX <- AFV5[!AFV5$f.924.0.0%in%c("Not walking"),]

# Fake the reference level "slow pace"
AFV5EX$f.924.0.0 <- relevel(AFV5EX$f.924.0.0, ref="Slow pace")
AFV5EX$f.924.0.0 <- factor(AFV5EX$f.924.0.0)


# Find subject included in all the imputation datasets #
IDcomT <- Reduce(intersect, list(AFT1EX$f.eid,AFT2EX$f.eid,AFT3EX$f.eid,AFT4EX$f.eid,AFT5EX$f.eid))
AFT1EX <- AFT1EX[AFT1EX$f.eid%in%IDcomT,]
AFT2EX <- AFT2EX[AFT2EX$f.eid%in%IDcomT,]
AFT3EX <- AFT3EX[AFT3EX$f.eid%in%IDcomT,]
AFT4EX <- AFT4EX[AFT4EX$f.eid%in%IDcomT,]
AFT5EX <- AFT5EX[AFT5EX$f.eid%in%IDcomT,]

IDcomV <- Reduce(intersect,list(AFV1EX$f.eid,AFV2EX$f.eid,AFV3EX$f.eid,AFV4EX$f.eid,AFV5EX$f.eid))
AFV1EX <- AFV1EX[AFV1EX$f.eid%in%IDcomV,]
AFV2EX <- AFV2EX[AFV2EX$f.eid%in%IDcomV,]
AFV3EX <- AFV3EX[AFV3EX$f.eid%in%IDcomV,]
AFV4EX <- AFV4EX[AFV4EX$f.eid%in%IDcomV,]
AFV5EX <- AFV5EX[AFV5EX$f.eid%in%IDcomV,]


#aa <- aggregate(All_riskFT,by=list(AFT1EX$age), mean)
#bb <- 1-FF095[26:56]

#pdf("female.pdf")
#plot(aa[,2],bb, xlab="5-year risk in UK biobank from prediction score #(log-scale)",ylab="5-year risk in UK from lifetables (log-scale)",xlim=c(0.005,0.09),ylim=c(0.005,0.09), pch="", log="xy")
#text(aa[,2],bb,labels=aa[,1], cex=0.8)
#abline(0,1)
#dev.off()

#aa <- aggregate(All_riskMT,by=list(AMT1EX$age), mean)
#bb <- 1-MM095[26:56]

#pdf("male.pdf")
#plot(aa[,2],bb, xlab="5-year risk in UK biobank from prediction score (log-scale)",ylab="5-year risk in UK from lifetables (log-scale)",xlim=c(0.008,0.12),ylim=c(0.008,0.12), pch="", log="xy")
#text(aa[,2],bb,labels=aa[,1], cex=0.8)
#abline(0,1)
#dev.off()


##### EXTRACT COEFFICIENTS AND FEANS ######

# Forumla
formufinF3 <- paste0("~",strsplit(as.character(MAINMODF$analyses[[1]]$formula),"~")[[3]])

# Coefficients
coxCC.temp1F <- coxph(as.formula(paste0("Surv(surv,out)",formufinF3)),data=AFT1EX)
coxCC.temp2F <- coxph(as.formula(paste0("Surv(surv,out)",formufinF3)),data=AFT2EX)
coxCC.temp3F <- coxph(as.formula(paste0("Surv(surv,out)",formufinF3)),data=AFT3EX)
coxCC.temp4F <- coxph(as.formula(paste0("Surv(surv,out)",formufinF3)),data=AFT4EX)
coxCC.temp5F <- coxph(as.formula(paste0("Surv(surv,out)",formufinF3)),data=AFT5EX)

ss1 <- summary(coxCC.temp1F)$coefficients[,1]
ss2 <- summary(coxCC.temp2F)$coefficients[,1]
ss3 <- summary(coxCC.temp3F)$coefficients[,1]
ss4 <- summary(coxCC.temp4F)$coefficients[,1]
ss5 <- summary(coxCC.temp5F)$coefficients[,1]

coefF <- rowMeans(cbind(ss1,ss2,ss3,ss4,ss5))

# Means
for (i in 1:length(list(VCFnewmean,GHFnewmean)))
	{coxCC.temp1F$means[grep(attributes(list(VCFnewmean,GHFnewmean)[[i]])$name,names(coxCC.temp1F$means))] <- list(VCFnewmean,GHFnewmean)[[i]][-1]}
mm1 <- coxCC.temp1F$means

for (i in 1:length(list(VCFnewmean,GHFnewmean)))
	{coxCC.temp2F$means[grep(attributes(list(VCFnewmean,GHFnewmean)[[i]])$name,names(coxCC.temp2F$means))] <- list(VCFnewmean,GHFnewmean)[[i]][-1]}
mm2 <- coxCC.temp1F$means

for (i in 1:length(list(VCFnewmean,GHFnewmean)))
	{coxCC.temp3F$means[grep(attributes(list(VCFnewmean,GHFnewmean)[[i]])$name,names(coxCC.temp3F$means))] <- list(VCFnewmean,GHFnewmean)[[i]][-1]}
mm3 <- coxCC.temp1F$means

for (i in 1:length(list(VCFnewmean,GHFnewmean)))
	{coxCC.temp4F$means[grep(attributes(list(VCFnewmean,GHFnewmean)[[i]])$name,names(coxCC.temp4F$means))] <- list(VCFnewmean,GHFnewmean)[[i]][-1]}
mm4 <- coxCC.temp1F$means

for (i in 1:length(list(VCFnewmean,GHFnewmean)))
	{coxCC.temp5F$means[grep(attributes(list(VCFnewmean,GHFnewmean)[[i]])$name,names(coxCC.temp5F$means))] <- list(VCFnewmean,GHFnewmean)[[i]][-1]}
mm5 <- coxCC.temp1F$means

meansF <- rowMeans(cbind(mm1,mm2,mm3,mm4,mm5))

# Variance-covariance matrix - same as pool(model)$ubar
vs1 <- vcov(coxCC.temp1F)
vs2 <- vcov(coxCC.temp2F)
vs3 <- vcov(coxCC.temp3F)
vs4 <- vcov(coxCC.temp4F)
vs5 <- vcov(coxCC.temp5F)

x <- array(c(vs1,vs2,vs3,vs4,vs5), dim=c(dim(vs1)[1],dim(vs1)[2],5)) 
vcovF <- apply(x,c(1,2),mean)

# Write data (ordered)
ordetobe <- c("age","f.2734.0.02","f.2734.0.03","f.2734.0.04","f.1239.0.02","f.1239.0.03","f.1249.0.02","f.1249.0.03","f.1249.0.04","f.1249.0.05","f.2178.0.02","f.2178.0.03","f.2178.0.04","f.2188.0.02","f.924.0.0Steady average pace","f.924.0.0Brisk pace","f.924.0.0None of the above","f.2090.0.02","f.2453.0.02","f.6145.0.l12","f.6145.0.l22","f.6145.0.l32","f.6145.0.l42","f.6145.0.l52","f.6145.0.l62","f.6145.0.l72","f.6146.0.l12","f.6146.0.l22","f.6146.0.l32","f.6146.0.l42","f.2453.0.02:age","age:f.6146.0.l12","age:f.6146.0.l32","age:f.6146.0.l42")

vcovFO <- vcovF[order(match(colnames(vs1),ordetobe)), order(match(colnames(vs1),ordetobe))]
coefFO <- coefF[order(match(names(coefF),ordetobe))]
meansFO <- meansF[order(match(names(meansF),ordetobe))]

# save 
write.csv(vcovFO,file="/proj/b2011036/uk.biobank/out_weights/varcovF.csv")
write.csv(cbind(coefFO,meansFO),file="/proj/b2011036/uk.biobank/out_weights/coefmeanF.csv")


###### BASELINE HAZARD ######
bb1 <- basehaz(coxCC.temp1F)$hazard[which(abs(basehaz(coxCC.temp1F)$time-5)==min(abs(basehaz(coxCC.temp1F)$time-5)))]
bb2 <-basehaz(coxCC.temp2F)$hazard[which(abs(basehaz(coxCC.temp2F)$time-5)==min(abs(basehaz(coxCC.temp2F)$time-5)))]
bb3 <-basehaz(coxCC.temp3F)$hazard[which(abs(basehaz(coxCC.temp3F)$time-5)==min(abs(basehaz(coxCC.temp3F)$time-5)))]
bb4 <-basehaz(coxCC.temp4F)$hazard[which(abs(basehaz(coxCC.temp4F)$time-5)==min(abs(basehaz(coxCC.temp4F)$time-5)))]
bb5 <-basehaz(coxCC.temp5F)$hazard[which(abs(basehaz(coxCC.temp5F)$time-5)==min(abs(basehaz(coxCC.temp5F)$time-5)))]

basehazF <- mean(c(bb1,bb2,bb3,bb4,bb5))

write.csv(basehazF,file="/proj/b2011036/uk.biobank/out_weights/basehazF.csv")


## Calculate the risk in the geographical validation ##
All_riskFT <- mclapply(seq(1:5), 		function(c) { calcweightedRISK(object=eval(parse(text=paste0("coxCC.temp",c,"F"))), newdata = eval(parse(text=paste0("AFV",c,"EX"))),newmeans=list(VCFnewmean,GHFnewmean),weights=wwFage)}, mc.cores=5)

All_riskF <- exp(-exp(rowMeans(cbind(log(-log(All_riskFT[[1]])),
log(-log(All_riskFT[[2]])),
log(-log(All_riskFT[[3]])),
log(-log(All_riskFT[[4]])),
log(-log(All_riskFT[[5]]))))))


save(All_riskF,file="/proj/b2011036/uk.biobank/prediction_results/All_riskF.Rdata")
save(AFV1EX,file="/proj/b2011036/uk.biobank/prediction_results/AFV1EX.Rdata")



### For supplementary table 5 ##

ind <- which(AFV1EX$sex=="Female" & AFV1EX$age==53 & AFV1EX$f.2734.0.0=="(1,2]" & AFV1EX$f.1239.0.0=="No" & AFV1EX$f.1249.0.0=="I have never smoked" & AFV1EX$f.2453.0.0=="No" & AFV1EX$f.2178.0.0=="Good" & AFV1EX$f.6145.0.l1=="None of the above" &  AFV1EX$f.924.0.0=="Slow pace" &  AFV1EX$f.6145.0.l1=="None of the above" & AFV1EX$f.2090.0.0=="No" & AFV1EX$f.2188.0.0=="No")
All_riskF[ind]
BIOAREF[ind,]



#################################################
#################################################
#### 8. BOOTSTRAP CI FOR CERTAIN QUANTITIES #####
#################################################
#################################################

k <- "f.2178.0.0"
CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"$surv"))), eval(parse(text=paste0("AMT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x)))), newdata=eval(parse(text=paste0("AMT",x)))))$concordance},mc.cores=5,mc.preschedule=FALSE))

CTbE <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"$surv"))), eval(parse(text=paste0("AMT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x)))), newdata=eval(parse(text=paste0("AMT",x)))))$std.err},mc.cores=5,mc.preschedule=FALSE))

mean(CTb)+1.96*mean(CTbE)
mean(CTb)-1.96*mean(CTbE)


k <- "f.2453.0.0"
CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"$surv"))), eval(parse(text=paste0("AFT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AFT",x)))), newdata=eval(parse(text=paste0("AFT",x)))))$concordance},mc.cores=5,mc.preschedule=FALSE))

CTbE <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"$surv"))), eval(parse(text=paste0("AFT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AFT",x)))), newdata=eval(parse(text=paste0("AFT",x)))))$std.err},mc.cores=5,mc.preschedule=FALSE))

mean(CTb)+1.96*mean(CTbE)
mean(CTb)-1.96*mean(CTbE)


k <- "f.924.0.0"
CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"$surv"))), eval(parse(text=paste0("AMT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x)))), newdata=eval(parse(text=paste0("AMT",x)))))$concordance},mc.cores=5,mc.preschedule=FALSE))

CTbE <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AMT",x,"$surv"))), eval(parse(text=paste0("AMT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age*",k)), eval(parse(text=paste0("AMT",x)))), newdata=eval(parse(text=paste0("AMT",x)))))$std.err},mc.cores=5,mc.preschedule=FALSE))

mean(CTb)+1.96*mean(CTbE)
mean(CTb)-1.96*mean(CTbE)


k <- "f.924.0.0"
CTb <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"$surv"))), eval(parse(text=paste0("AFT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AFT",x)))), newdata=eval(parse(text=paste0("AFT",x)))))$concordance},mc.cores=5,mc.preschedule=FALSE))

CTbE <- do.call("cbind",mclapply(1:5,function(x) {survConcordance(Surv(eval(parse(text=paste0("AFT",x,"$surv"))), eval(parse(text=paste0("AFT",x,"$out"))))~predict(coxph(as.formula(paste0("Surv(surv,out)~age+",k)), eval(parse(text=paste0("AFT",x)))), newdata=eval(parse(text=paste0("AFT",x)))))$std.err},mc.cores=5,mc.preschedule=FALSE))

mean(CTb)+1.96*mean(CTbE)
mean(CTb)-1.96*mean(CTbE)


### Likelihood ratio test for a model with only age and age+prediction model ##
class.temp1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT1)
class.temp2 <- coxph(Surv(surv,out==1)~age, AFT1)
anova(class.temp2,class.temp1)	

### Likelihood ratio test for a model with charlson index + prediction score ##
class.temp1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3)), AFT1)
class.temp2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF3,"+as.numeric(charlsonSR)")), AFT1)
anova(class.temp2,class.temp1)	


### BOOTSTRAP FOR DIFFERENCE BETWEEN AGE + CHARSON INDEX AND AGE + PREDICTIN SCORE ####
### WOMEN ###

formufinF <- paste0("~",strsplit(as.character(MAINMODF$analyses[[1]]$formula),"~")[[3]])
formufinF5 <- paste0("~","age+charlsonSR")
m1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT1)
m2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT1)

CC <- NULL
for(i in 1:200)
{
	set.seed(123+i)
	ind <- sample(1:nrow(AFV1), replace=T)
	
	cc <- survConcordance(Surv(AFV1[ind,]$surv,AFV1[ind,]$out)~predict(m1, AFV1[ind,]))$concordance - survConcordance(Surv(AFV1[ind,]$surv,AFV1[ind,]$out)~predict(m2, AFV1[ind,]))$concordance
	print(i)
  CC <- c(CC,cc)	
}

2*pnorm(-abs(mean(CC)/sd(CC)))



### BOOTSTRAP FOR DIFFERENCE BETWEEN AGE + CHARSON INDEX AND AGE + PREDICTIN SCORE ####
#### MEN ####


formufinM <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]])
formufinM5 <- paste0("~","age+charlsonSR")
m1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT1)
m2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT1)

CC <- NULL
for(i in 1:200)
{
	set.seed(123+i)
	ind <- sample(1:nrow(AMV1), replace=T)
	
	cc <- survConcordance(Surv(AMV1[ind,]$surv,AMV1[ind,]$out)~predict(m1, AMV1[ind,]))$concordance - survConcordance(Surv(AMV1[ind,]$surv,AMV1[ind,]$out)~predict(m2, AMV1[ind,]))$concordance
	print(i)
  CC <- c(CC,cc)	
}

2*pnorm(-abs(mean(CC)/sd(CC)))




########################################
########################################
#### 9. EXTRA ANALYSES FOR REVIWERS ####
########################################
########################################

##########################################################################################################
## 1. CHECK THIS RESULTS IN COMPARISON WITH OTHER ANOTHER STUDY WHICH IS FROM SURVEY DATA PMID: 9883792 ##
##########################################################################################################

AM1$f.2178.0.0_dic <- ifelse(AM1$f.2178.0.0%in%c("Poor","Fair"),1,0)
AM2$f.2178.0.0_dic <- ifelse(AM2$f.2178.0.0%in%c("Poor","Fair"),1,0)
AM3$f.2178.0.0_dic <- ifelse(AM3$f.2178.0.0%in%c("Poor","Fair"),1,0)
AM4$f.2178.0.0_dic <- ifelse(AM4$f.2178.0.0%in%c("Poor","Fair"),1,0)
AM5$f.2178.0.0_dic <- ifelse(AM5$f.2178.0.0%in%c("Poor","Fair"),1,0)


modb <- tryCatch(
	coxph.mids(as.formula(paste0("Surv(age,age_exit,out)~f.2178.0.0_dic")), 5, N_MISSM,"AM"),
	error=function(e) {NULL})

pp <- pool(modb)
cbind(levels(AM1[,k]),c("Reference",pp$qbar),c("Reference",sqrt(diag(pp$t))))

exp(pp$qbar)
exp(pp$qbar + 1.96*sqrt(diag(pp$t)))
exp(pp$qbar - 1.96*sqrt(diag(pp$t)))
levels(AM1$f.2178.0.0_dic)

############################################
### 2. FRAILTY MODELS - Reviewer Table 1####
############################################

### MEN ####
load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
AM1$center2 <- factor(AM1$center)

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")


RESFM <- NULL
RESBM <- NULL
for (i in all.vars(MAINMODM$analyses[[1]]$formula)[c(-1,-2,-4)])
{
	tt <- table(AM1[,i])
	AM1[,i]  <- relevel(AM1[,i], ref=names(tt)[tt==max(tt)][1])
	modF <- coxph(Surv(age,age_exit,out)~AM1[,i] + frailty.gamma(center2), data=AM1)
	modb <- coxph(Surv(age,age_exit,out)~AM1[,i], data=AM1)	
	resF <- paste0(formatC(exp(coef(modF)),digits=1,format="f")," [",formatC(exp(coef(modF)-1.96*coef(summary(modF))[1:nrow(coef(summary(modF)))-1,2]), digits=1, format="f"),"-",formatC(exp(coef(modF)+1.96*coef(summary(modF))[1:nrow(coef(summary(modF)))-1,2]),digits=1, format="f"),"]")
	resB <- paste0(formatC(exp(coef(modb)),digits=1,format="f")," [",formatC(exp(coef(modb)-1.96*coef(summary(modb))[,3]),digits=1,format="f"),"-",formatC(exp(coef(modb)+1.96*coef(summary(modb))[,3]),digits=1,format="f"),"]")
	RESFM <- rbind(RESFM,cbind(i,resF))
	RESBM <- rbind(RESBM,cbind(i,resB))	
	print(i)
}

write.csv(cbind(RESFM,RESBM),file="test.csv")



#############################################
#### 3. RANKING-BASED OUTCOME COMPARISON ####
#############################################


## MEN ##
resM <- read.xls("/proj/b2011036/uk.biobank/univariateM.xlsx", header=T, stringsAsFactor=F)

### PROCESS DATA IN MEN ###
for (dis in c("all","Healthy","CA","CVD","RE","DG","EC","OD"))
{
	K <- NULL
	for (f in unique(resM$Original.code))
	{
		
		t <- resM[resM$Original.code == f,]
		if(unique(t$P.value.for.Schoenfeld.residuals)!= "&lt;0.00001")
		{
			t2 <- t[t$Age.category=="All",]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- unique(eval(parse(text=paste0("t2$C.index.",dis))))
					}	
					else {k <- NA}
							
		}	
		else if(unique(t$P.value.for.Schoenfeld.residuals)== "&lt;0.00001")
		{
			t2 <- t[!t$Age.category%in%c("All","Reference"),]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- unique(eval(parse(text=paste0("t2$C.index.",dis))))
					}		
					else {k <- NA}		
		}
		
		K <- rbind(K,c(f,k,unique(t$C4),unique(t$UK.Biobank.code)))
	}	
	
	assign(paste0("MCIN",dis),K)	
}


# Merge the data
a1 <- merge(data.frame(C.index.all=as.numeric(MCINall[,2]), Code=MCINall[,1],UK.Biobank.code=MCINall[,4] , Measurement=MCINall[,3] ),data.frame(C.index.CA=as.numeric(MCINCA[,2]),Code=MCINall[,1]), by="Code")
a2 <- merge(a1,data.frame(C.index.CVD=as.numeric(MCINCVD[,2]), Code=MCINCVD[,1]), by="Code")
a3 <- merge(a2,data.frame(C.index.RE=as.numeric(MCINRE[,2]), Code=MCINRE[,1]), by="Code")
a4 <- merge(a3,data.frame(C.index.DG=as.numeric(MCINDG[,2]), Code=MCINDG[,1]), by="Code")
a5 <- merge(a4,data.frame(C.index.EC=as.numeric(MCINEC[,2]), Code=MCINEC[,1]), by="Code")
a6 <- merge(a5,data.frame(C.index.OD=as.numeric(MCINOD[,2]), Code=MCINOD[,1]), by="Code")

# keep only variables without missing or lack of convergence
RANKM <- a6[rowSums(is.na(a6[,5:ncol(a6)]))==0,]
RANKM2 <- apply(-RANKM[,5:ncol(RANKM)],2,rank)

## Global ranking ##
global <- rank(rowSums(RANKM2), ties.method="random")
global_sd <- apply(RANKM2,1,sd)
global_annM <- cbind(as.character(RANKM[,3]),as.character(RANKM[,4]), rank(-RANKM[,2], ties.method="random") ,global,global_sd)



## WOMEN ##

resF <- read.xls("/proj/b2011036/uk.biobank/univariateF.xlsx", header=T, stringsAsFactor=F)

### PROCESS DATA IN MEN ###
for (dis in c("all","Healthy","CA","CVD","RE","DG","EC","OD"))
{
	K <- NULL
	for (f in unique(resF$Original.code))
	{
		
		t <- resF[resF$Original.code == f,]
		if(unique(t$P.value.for.Schoenfeld.residuals)!= "&lt;0.00001")
		{
			t2 <- t[t$Age.category=="All",]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- unique(eval(parse(text=paste0("t2$C.index.",dis))))
					}	
					else {k <- NA}
							
		}	
		else if(unique(t$P.value.for.Schoenfeld.residuals)== "&lt;0.00001")
		{
			t2 <- t[!t$Age.category%in%c("All","Reference"),]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- unique(eval(parse(text=paste0("t2$C.index.",dis))))
					}		
					else {k <- NA}		
		}
		
		K <- rbind(K,c(f,k,unique(t$C4),unique(t$UK.Biobank.code)))
	}	
	
	assign(paste0("FCIN",dis),K)	
}


# Merge the data
a1 <- merge(data.frame(C.index.all=as.numeric(FCINall[,2]), Code=FCINall[,1],UK.Biobank.code=FCINall[,4] , Measurement=FCINall[,3] ),data.frame(C.index.CA=as.numeric(FCINCA[,2]),Code=FCINall[,1]), by="Code")
a2 <- merge(a1,data.frame(C.index.CVD=as.numeric(FCINCVD[,2]), Code=FCINCVD[,1]), by="Code")
a3 <- merge(a2,data.frame(C.index.RE=as.numeric(FCINRE[,2]), Code=FCINRE[,1]), by="Code")
a4 <- merge(a3,data.frame(C.index.DG=as.numeric(FCINDG[,2]), Code=FCINDG[,1]), by="Code")
a5 <- merge(a4,data.frame(C.index.EC=as.numeric(FCINEC[,2]), Code=FCINEC[,1]), by="Code")
a6 <- merge(a5,data.frame(C.index.OD=as.numeric(FCINOD[,2]), Code=FCINOD[,1]), by="Code")

# keep only variables without missing or lack of convergence
RANKF <- a6[rowSums(is.na(a6[,5:ncol(a6)]))==0,]
RANKF2 <- apply(-RANKF[,5:ncol(RANKF)],2,rank)

## Global ranking ##
global <- rank(rowSums(RANKF2), ties.method="random")
global_sd <- apply(RANKF2,1,sd)
global_annF <- cbind(as.character(RANKF[,3]),as.character(RANKF[,4]), rank(-RANKF[,2], ties.method="random") ,global,global_sd)

## MERGE MEN AND WOMEN ##
global_annMF <- merge(data.frame(global_annM),data.frame(global_annF), by="V1", all.x=T, all.y=T)
write.csv(global_annMF,file="/proj/b2011036/uk.biobank/out_weights/ranking_for_suppl_table3.csv")


################################################
### 4. EXCLUDE SOCIO-DEMOGRAPHICS INDICATORS ###
################################################

###########
### MEN ###
###########

NAMESM <- colnames(AM1)[!colnames(AM1)%in%c("age","age_exit","sex","surv","out","f.eid","center","nal","surv5","out5","status","cofd","agecat","charlsonSR")]

resM <- read.xls("/proj/b2011036/uk.biobank/univariateM.xlsx", header=T, stringsAsFactor=F)

##EXCLUDE ALL CATEGORIES SHOULD NOT BE INCLUDED IN THE SCORE ##
resM2 <- resM[!resM$Measurement.class%in%c("Physical measures","Blood assays","Cognitive function","Sociodemographics"),]


## EXTRACT TOP 20 MEASUREMENTS FOR EACH CAUSE OF DEATH
FFM20 <- NULL
for (dis in c("CA","CVD","RE","DG","EC","OD"))
{
	K <- NULL
	for (f in unique(resM2$Original.code))
	{
		
		t <- resM2[resM2$Original.code == f,]
		if(unique(t$P.value.for.Schoenfeld.residuals)!= "&lt;0.00001")
		{
			t2 <- t[t$Age.category=="All",]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),unique(t$Original.code),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}	
					else {k <- c(unique(t$Original.code),unique(t$Original.code),NA)}
							
		}	
		else if(unique(t$P.value.for.Schoenfeld.residuals)== "&lt;0.00001")
		{
			t2 <- t[!t$Age.category%in%c("All","Reference"),]
					if (unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "Singular fit" & unique(eval(parse(text=paste0("t2$C.index.",dis)))) != "N &lt;= 5")
					{
						k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),unique(eval(parse(text=paste0("t2$C.index.",dis)))))
					}		
					else {k <- c(unique(t$Original.code),paste0(unique(t$Original.code),"*age"),NA)}		
		}
		
		K <- rbind(K,k)
	}	
	
	FFM20 <- cbind(FFM20,K[order(-as.numeric(K[,3])),2][1:20])
}



to_sel <- unique(as.vector(FFM20))
formufinM <- paste0("~",paste0(to_sel, collapse="+"))

### VARIABLE SELECTION ###

# These should be excluded because from verbal interview
to_ex <- c("f.134.0.0*age","f.135.0.0*age","f.137.0.0*age","f.189.0.0","f.6162.0.l4*age")
to_selM <- to_sel[!to_sel%in%to_ex]
formufinM2 <- paste0("~",paste0(to_selM, collapse="+"))

set.seed(123)
MODSTEPW <- cph(as.formula(paste0("Surv(surv,out==1)",formufinM2)),data=AMT1,x=T,y=T)
STEPWM <- fastbw(MODSTEPW,type="residual")

ph <- read.xls("/proj/b2011036/uk.biobank/coded_questions.xlsx", header=T, stringsAsFactor=F)
ph[as.character(ph[,1])%in%unlist(lapply(strsplit(STEPWM$names.kept,"\\."),"[",2)),c(6,7)]

# Add extra variables for conformity with uk biobank # 
newaddM <- STEPWM$names.kept

MAINMODMR <- coxph.mids(as.formula(paste0("Surv(surv,out==1)~",paste0(newaddM,collapse="+"))), 5, N_MISSM,"AMT")



# DISCRIMINATION #
formufinMR <- paste0("~",strsplit(as.character(MAINMODMR$analyses[[1]]$formula),"~")[[3]])

class.temp1R <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT1)			
class.temp2R <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT2)
class.temp3R <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT3)
class.temp4R <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT4)
class.temp5R <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT5)

c1R <- survConcordance(Surv(AMV1$surv,AMV1$out==1)~predict(class.temp1R,AMV1))$concordance
c2R <- survConcordance(Surv(AMV2$surv,AMV2$out==1)~predict(class.temp2R,AMV2))$concordance
c3R <- survConcordance(Surv(AMV3$surv,AMV3$out==1)~predict(class.temp3R,AMV3))$concordance
c4R <- survConcordance(Surv(AMV4$surv,AMV4$out==1)~predict(class.temp4R,AMV4))$concordance
c5R <- survConcordance(Surv(AMV5$surv,AMV5$out==1)~predict(class.temp5R,AMV5))$concordance

CINDEXMR <- mean(c(c1R,c2R,c3R,c4R,c5R))


### Bootstrapp for differences ##
formufinM <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]])

m1 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinMR)), AMT1)
m2 <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT1)

CC <- NULL
for(i in 1:200)
{
	set.seed(123+i)
	ind <- sample(1:nrow(AMV1), replace=T)
	
	cc <- survConcordance(Surv(AMV1[ind,]$surv,AMV1[ind,]$out)~predict(m1, AMV1[ind,]))$concordance - survConcordance(Surv(AMV1[ind,]$surv,AMV1[ind,]$out)~predict(m2, AMV1[ind,]))$concordance
	print(i)
  CC <- c(CC,cc)	
}

2*pnorm(-abs(mean(CC)/sd(CC)))



###############################################################################
### 5. COMPARISON WITH CHARLSON SCORE IN INDIVIDUAL WITH CHARLSON SCORE > 0 ###
###############################################################################

#### MEN ###
## Prediction model #
load("/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")
formufinM <- paste0("~",strsplit(as.character(MAINMODM$analyses[[1]]$formula),"~")[[3]])

class.temp1CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT1[AMT1$charlsonSR != 0,])			
class.temp2CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT2[AMT2$charlsonSR != 0,])
class.temp3CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT3[AMT3$charlsonSR != 0,])
class.temp4CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT4[AMT4$charlsonSR != 0,])
class.temp5CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM)), AMT5[AMT5$charlsonSR != 0,])

c1CC <- survConcordance(Surv(AMV1[AMV1$charlsonSR != 0,]$surv,AMV1[AMV1$charlsonSR != 0,]$out==1)~predict(class.temp1CC,AMV1[AMV1$charlsonSR != 0,]))$concordance
c2CC <- survConcordance(Surv(AMV2[AMV2$charlsonSR != 0,]$surv,AMV2[AMV2$charlsonSR != 0,]$out==1)~predict(class.temp2CC,AMV2[AMV2$charlsonSR != 0,]))$concordance
c3CC <- survConcordance(Surv(AMV3[AMV3$charlsonSR != 0,]$surv,AMV3[AMV3$charlsonSR != 0,]$out==1)~predict(class.temp3CC,AMV3[AMV3$charlsonSR != 0,]))$concordance
c4CC <- survConcordance(Surv(AMV4[AMV4$charlsonSR != 0,]$surv,AMV4[AMV4$charlsonSR != 0,]$out==1)~predict(class.temp4CC,AMV4[AMV4$charlsonSR != 0,]))$concordance
c5CC <- survConcordance(Surv(AMV5[AMV5$charlsonSR != 0,]$surv,AMV5[AMV5$charlsonSR != 0,]$out==1)~predict(class.temp5CC,AMV5[AMV5$charlsonSR != 0,]))$concordance

CINDEXCCM <- mean(c(c1CC,c2CC,c3CC,c4CC,c5CC))

## Age + Charlson index ##

formufinM5 <- paste0("~","age+charlsonSR")

## Age + Charlson score ##
class.temp1CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT1[AMT1$charlsonSR != 0,])			
class.temp2CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT2[AMT2$charlsonSR != 0,])
class.temp3CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT3[AMT3$charlsonSR != 0,])
class.temp4CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT4[AMT4$charlsonSR != 0,])
class.temp5CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinM5)), AMT5[AMT5$charlsonSR != 0,])

c1CC <- survConcordance(Surv(AMV1[AMV1$charlsonSR != 0,]$surv,AMV1[AMV1$charlsonSR != 0,]$out==1)~predict(class.temp1CC,AMV1[AMV1$charlsonSR != 0,]))$concordance
c2CC <- survConcordance(Surv(AMV2[AMV2$charlsonSR != 0,]$surv,AMV2[AMV2$charlsonSR != 0,]$out==1)~predict(class.temp2CC,AMV2[AMV2$charlsonSR != 0,]))$concordance
c3CC <- survConcordance(Surv(AMV3[AMV3$charlsonSR != 0,]$surv,AMV3[AMV3$charlsonSR != 0,]$out==1)~predict(class.temp3CC,AMV3[AMV3$charlsonSR != 0,]))$concordance
c4CC <- survConcordance(Surv(AMV4[AMV4$charlsonSR != 0,]$surv,AMV4[AMV4$charlsonSR != 0,]$out==1)~predict(class.temp4CC,AMV4[AMV4$charlsonSR != 0,]))$concordance
c5CC <- survConcordance(Surv(AMV5[AMV5$charlsonSR != 0,]$surv,AMV5[AMV5$charlsonSR != 0,]$out==1)~predict(class.temp5CC,AMV5[AMV5$charlsonSR != 0,]))$concordance

CINDEXCCCM <- mean(c(c1CC,c2CC,c3CC,c4CC,c5CC))



#### WOMEN ###
## Prediction model #
load("/proj/b2011036/uk.biobank/prediction_results/MAINMODF.Rdata")
formufinF <- paste0("~",strsplit(as.character(MAINMODF$analyses[[1]]$formula),"~")[[3]])

class.temp1CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT1[AFT1$charlsonSR != 0,])			
class.temp2CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT2[AFT2$charlsonSR != 0,])
class.temp3CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT3[AFT3$charlsonSR != 0,])
class.temp4CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT4[AFT4$charlsonSR != 0,])
class.temp5CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF)), AFT5[AFT5$charlsonSR != 0,])

c1CC <- survConcordance(Surv(AFV1[AFV1$charlsonSR != 0,]$surv,AFV1[AFV1$charlsonSR != 0,]$out==1)~predict(class.temp1CC,AFV1[AFV1$charlsonSR != 0,]))$concordance
c2CC <- survConcordance(Surv(AFV2[AFV2$charlsonSR != 0,]$surv,AFV2[AFV2$charlsonSR != 0,]$out==1)~predict(class.temp2CC,AFV2[AFV2$charlsonSR != 0,]))$concordance
c3CC <- survConcordance(Surv(AFV3[AFV3$charlsonSR != 0,]$surv,AFV3[AFV3$charlsonSR != 0,]$out==1)~predict(class.temp3CC,AFV3[AFV3$charlsonSR != 0,]))$concordance
c4CC <- survConcordance(Surv(AFV4[AFV4$charlsonSR != 0,]$surv,AFV4[AFV4$charlsonSR != 0,]$out==1)~predict(class.temp4CC,AFV4[AFV4$charlsonSR != 0,]))$concordance
c5CC <- survConcordance(Surv(AFV5[AFV5$charlsonSR != 0,]$surv,AFV5[AFV5$charlsonSR != 0,]$out==1)~predict(class.temp5CC,AFV5[AFV5$charlsonSR != 0,]))$concordance

CINDEXCCF <- mean(c(c1CC,c2CC,c3CC,c4CC,c5CC))

## Age + Charlson index ##

formufinF5 <- paste0("~","age+charlsonSR")

## Age + Charlson score ##
class.temp1CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT1[AFT1$charlsonSR != 0,])			
class.temp2CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT2[AFT2$charlsonSR != 0,])
class.temp3CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT3[AFT3$charlsonSR != 0,])
class.temp4CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT4[AFT4$charlsonSR != 0,])
class.temp5CC <- coxph(as.formula(paste0("Surv(surv,out==1)",formufinF5)), AFT5[AFT5$charlsonSR != 0,])

c1CC <- survConcordance(Surv(AFV1[AFV1$charlsonSR != 0,]$surv,AFV1[AFV1$charlsonSR != 0,]$out==1)~predict(class.temp1CC,AFV1[AFV1$charlsonSR != 0,]))$concordance
c2CC <- survConcordance(Surv(AFV2[AFV2$charlsonSR != 0,]$surv,AFV2[AFV2$charlsonSR != 0,]$out==1)~predict(class.temp2CC,AFV2[AFV2$charlsonSR != 0,]))$concordance
c3CC <- survConcordance(Surv(AFV3[AFV3$charlsonSR != 0,]$surv,AFV3[AFV3$charlsonSR != 0,]$out==1)~predict(class.temp3CC,AFV3[AFV3$charlsonSR != 0,]))$concordance
c4CC <- survConcordance(Surv(AFV4[AFV4$charlsonSR != 0,]$surv,AFV4[AFV4$charlsonSR != 0,]$out==1)~predict(class.temp4CC,AFV4[AFV4$charlsonSR != 0,]))$concordance
c5CC <- survConcordance(Surv(AFV5[AFV5$charlsonSR != 0,]$surv,AFV5[AFV5$charlsonSR != 0,]$out==1)~predict(class.temp5CC,AFV5[AFV5$charlsonSR != 0,]))$concordance

CINDEXCCCF <- mean(c(c1CC,c2CC,c3CC,c4CC,c5CC))


######################################
### % OF MISSING FOR TOP VARIABLES ###
######################################

load("/proj/b2011036/uk.biobank/prediction_results/MAINMODM.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/MAINMODF.Rdata")
load("/proj/b2011036/uk.biobank/out5.Rdata")

subj_miss <- rowSums(is.na(bdE5))
bdE6 <- bdE5[(subj_miss/ncol(bdE5)) < 0.8,]

# Get top 20 variables and those included in the prediction model #

VARMF <- unique(c(NAMESM[FFM%in%FFM20], NAMESF[FFF%in%FFF20], all.vars(MAINMODM$analyses[[1]]$formula),all.vars(MAINMODF$analyses[[1]]$formula)))[c(-87,-88,-89)]

## PLOT N. MISSING PER VARIABLE ##
var_miss <- colSums(is.na(bdE6))
sort(var_miss[colnames(bdE6)%in%VARMF]/nrow(bdE6))
summary(sort(var_miss[colnames(bdE6)%in%VARMF]/nrow(bdE6))[1:97])

