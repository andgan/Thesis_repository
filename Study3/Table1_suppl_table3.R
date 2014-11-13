#----------------------------------------------
# Filename: Table1_suppl_table3.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 11JUL2014 - Add LDL-cholesterol and triglycerides as covariates. Add a analysis not included in the paper that studies these associations in PIVUS
#          23SEP2014 - Added analysis requested by the reviewer - Fatal vs. non fatal CHD.
# Purpose: Analysis used in the paper. Table 1 (meta-analysis results) and supplementary Table 3 (Single-study results)
# Note: Association between top-findings and incidenct CHD, adjusted for FHS.
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small)
# Data created: 
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/


### FUNCTION
X.coxph <- 
function (object, newdata, na.action = na.pass)
{
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      #.checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    X <- X[,-1, drop=FALSE]
  }
  X - rep(object$means, each = nrow(X)) # !!!
}



#### LOAD ULSAM ####

library(survival)

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/ulsam.metabolomics/data_processed/plasma/ulsam_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]


metabo_sub <- metabo_p[metabo_p$time==0  & metabo_p$double_==0,]


##### SELECT CHD EVENTS #####
metabo_sub_nopchdU <- metabo_sub[metabo_sub$chd==0 | (metabo_sub$chd==1 & metabo_sub$chdd > metabo_sub$check_d),]

#### 10 years follow-up #####
metabo_sub_nopchdU$incchd10 <- ifelse(metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd<=10 & metabo_sub_nopchdU$incchd==1,1,0)
metabo_sub_nopchdU$surv_chd10 <- metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd
metabo_sub_nopchdU$surv_chd10[metabo_sub_nopchdU$surv_chd10>10] <- 10
metabo_sub_nopchdU$surv_chd10[metabo_sub_nopchdU$surv_chd10<0] <- 0.0001
metabo_sub_nopchdU$surv_chd <- metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd

## Log-crp
metabo_sub_nopchdU$logcrp <- log(metabo_sub_nopchdU$crp)
metabo_sub_nopchdU$logtg <- log(metabo_sub_nopchdU$tg)


#### 20 years follow-up #####
metabo_sub_nopchdU$incchd20 <- ifelse(metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd<=20 & metabo_sub_nopchdU$incchd==1,1,0)
metabo_sub_nopchdU$surv_chd20 <- metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd
metabo_sub_nopchdU$surv_chd20[metabo_sub_nopchdU$surv_chd20>20] <- 20
metabo_sub_nopchdU$surv_chd20[metabo_sub_nopchdU$surv_chd20<0] <- 0.0001


cor(metabo_sub$glucose,as.numeric(metabo_sub$M225.035T32.769), use="complete.obs")


##### LOAD TWINGENE #####

load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/twge.metabolomics/data_processed/twge_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]

### Recode diabetes ###
metabo_p$diab_baseline_new <- ifelse((metabo_p$diabetess==1 & metabo_p$diabetess_doctor==1) | metabo_p$diabetess_medicine==1 | (metabo_p$glucos>= 7 & metabo_p$fasting==0) | (metabo_p$glucos>= 11 & metabo_p$fasting!=0) | (metabo_p$diab_regd <= metabo_p$check_d & !is.na(metabo_p$diab_regd)),1,0)

metabo_sub_chd <- metabo_p[(metabo_p$sub==1  | metabo_p$incchd==1) & metabo_p$fasting!=2,]



# CHD-free subjects
metabo_sub_nopchdT <- metabo_sub_chd[metabo_sub_chd$chd==0 | (metabo_sub_chd$chd==1 & metabo_sub_chd$chdd > metabo_sub_chd$check_d),]

## Log-crp
metabo_sub_nopchdT$logcrp <- log(metabo_sub_nopchdT$crp)
metabo_sub_nopchdT$logtg <- log(metabo_sub_nopchdT$tg)


library(Hmisc)


### ULSAM ###

## TF MG 18:2, TF monosaccarides, TF LysoPC 18:2. 
## 

topfU <- c("M393.238T397.936","M225.035T32.769","M520.340T357.723","M221.081T145.997","M661.528T590.990","M552.403T462.976","M570.356T377.717","M542.325T331.776")



RESU <- NULL
RESU2 <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU[,i]))+sbp+bmi+smoke+antihyp+ldl+hdl+diab+log(tg), data =metabo_sub_nopchdU)
	res <- c(i,
		paste(round(exp(as.numeric(mod$coefficients[2])),2)," (",
		round(exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),2),"-",
		round(exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]),2),")",sep=""), 	2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])	
	RESU <- rbind(RESU,res)
	RESU2 <- rbind(RESU2,res2)
}





#### TWINGENE ALL ####

topfT <- c("M393.231T384.578","M225.035T32.546","M520.340T346.671","M221.081T142.576","M661.528T580.694","M552.403T449.740",
"M570.356T366.208","M542.325T321.861")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$tg) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi),]


metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1


## Age in 3 categories ###
metabo_sub_nopchdT4$age3 <- cut(metabo_sub_nopchdT4$age, breaks=c(min(metabo_sub_nopchdT4$age),65,70,max(metabo_sub_nopchdT4$age)),include.lowest = T)

STRATA <- NULL
STRATAM <- NULL
INT <- NULL
MAIN <- NULL
MAIN2 <- NULL
for (i in topfT)
{

	metabo_sub_nopchdT4$x <- scale(as.numeric(metabo_sub_nopchdT4[,i]))
	
	#### CHECK INTERACTION WITH AGE ###
	stratcox2<-coxph(Surv(surv_chd,incchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox2,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox2$var)[1],dim(stratcox2$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox2$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[12])/as.numeric(adjse[12])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(surv_chd,incchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

	dfb<-as.matrix(resid(stratcox3,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox3$var)[1],dim(stratcox3$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox3$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	MAIN <- rbind(MAIN,c(i,paste(round(exp(as.numeric(stratcox3$coefficients[3])),2)," (",round(exp(stratcox3$coefficients[3]+qnorm(0.025)*adjse[3]),2),"-",round(exp(stratcox3$coefficients[3]-qnorm(0.025)*adjse[3]),2),")",sep=""), 2*pnorm(-abs(as.numeric(stratcox3$coefficients[3])/as.numeric(adjse[3])))))
 MAIN2 <- rbind(MAIN2,c(i,as.numeric(stratcox3$coefficients[3]),adjse[3]))


 #### STRATIFIED ####
	stratcox<-coxph(Surv(surv_chd,incchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox$var+gamma


 ###  LOOP OVER DIFFERENT AGE-LEVELS ####
	RES1 <- NULL
	RES2 <- NULL
  for (k in 1:length(levels(metabo_sub_nopchdT4$age3)))
	 {
		newdata0 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=0,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		RES1 <- rbind(RES1,res1)
		res2 <- c(i,lp,se)
		RES2 <- rbind(RES2,res2)
	 }
	
	# Report stratified analysis
	STRATA <- rbind(STRATA,RES1)
	STRATAM <- rbind(STRATAM,RES2[3,])
	
}

## Check if significant interaction
INT
## Check MAIN effect results
MAIN
## Check MAIN effect for meta-analysis
MAIN2
## Check resutls for stratified analysis
STRATA
### Check results for stratified analysis (age > 70), result for meta-analysis
STRATAM


#### Only if interaction use the STRATAM, otherwise just use MAIN ####
#### SO use STRATAM only for M520.340T346.671 ####
REST2 <- rbind(MAIN2[1,],MAIN2[2,],STRATAM[3,],MAIN2[4,],MAIN2[5,],MAIN2[6,],MAIN2[7,],MAIN2[8,])




############################
###### META-ANALYSIS #######
############################

library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:length(REST2[,1]))
{
	beta <- c(as.numeric(REST2[i,2]),as.numeric(RESU2[i,2]))
	se <- c(as.numeric(REST2[i,3]),as.numeric(RESU2[i,3]))
	temp <- metagen(beta,se,studlab=c("TwinGene","Ulsam"))
	meta_b <- c(meta_b,temp$TE.random)
	meta_s <- c(meta_s,temp$seTE.random)
}

chd_meta <- cbind(topfT,paste(round(exp(meta_b),2)," (",round(exp(meta_b+qnorm(0.025)*meta_s),2),"-",round(exp(meta_b-qnorm(0.025)*meta_s),2),")",sep=""),2*pnorm(-abs(meta_b/meta_s)))





###################################################
#### REVIEWER ANALYSIS: FATAL Vs. NON-FATAL MI ####
###################################################


#############
### FATAL ###
#############


##### ULSAM #####

library(survival)

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/ulsam.metabolomics/data_processed/plasma/ulsam_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]


metabo_sub <- metabo_p[metabo_p$time==0  & metabo_p$double_==0,]

## Now read mortality data ##
library(foreign) 
dd <- read.dta("/home/andrea/glob/metabo_chd/death111013.dta")
dd$p013 <- as.Date(dd$p013,"%Y%m%d")
dd <- dd[,c("pat","chd_death","p013","p002")]

# Merge with mortality data
metabo_sub2 <- merge(metabo_sub,dd,by="pat")
# If missing then assign end of follow-up
metabo_sub2$p013[is.na(metabo_sub2$p013)] <- as.Date("20110831","%Y%m%d")
# Survival
metabo_sub2$death_surv <- (metabo_sub2$p013-metabo_sub2$check_d)/365.25
# Truncate at 10 years
metabo_sub2$death_surv[metabo_sub2$death_surv>10] <- 10
metabo_sub2$chd_death <- ifelse(metabo_sub2$death_surv<=10 & metabo_sub2$chd_death==1,1,0)

# Option 1
# Exclude individuals with non-fatal CHD event #
metabo_sub_nopchdU <- metabo_sub2[metabo_sub2$chd==0 | (metabo_sub2$chd_death==1 & metabo_sub2$incchd==1 & metabo_sub2$p013==metabo_sub2$chdd & (metabo_sub2$p013-metabo_sub2$check_d)/365.25 <=10),]




# Indicator if the mortality the incident CHD is due to fatal MI
twexcl <- 	metabo_sub2$pat[metabo_sub2$chd_death==1 & metabo_sub2$p013==metabo_sub2$chdd & !is.na(metabo_sub2$p013) & (metabo_sub2$p013-metabo_sub2$check_d)/365.25 <=10]

metabo_sub_nopchdU <- metabo_sub2[(metabo_sub2$chd==0 | (metabo_sub2$chd==1 & metabo_sub2$chdd > metabo_sub2$check_d)) & !metabo_sub2$pat%in%twexcl,]





# Option 2
# Keep individuals with non-fatal CHD event #
metabo_sub_nopchdU <- metabo_sub2


topfU <- c("M393.238T397.936","M520.340T357.723","M661.528T590.990","M522.356T395.236")

RESU <- NULL
RESU2 <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(as.numeric(death_surv),chd_death) ~ age+scale(as.numeric(metabo_sub2[,i]))+sbp+bmi+smoke+antihyp+ldl+hdl+diab+log(tg), data =metabo_sub2)
	res <- c(i,
		paste(round(exp(as.numeric(mod$coefficients[2])),2)," (",
		round(exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),2),"-",
		round(exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]),2),")",sep=""), 	2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])	
	RESU <- rbind(RESU,res)
	RESU2 <- rbind(RESU2,res2)
}



##### TWINGENE #####

library(survival)

load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/twge.metabolomics/data_processed/twge_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]

### Recode diabetes ###
metabo_p$diab_baseline_new <- ifelse((metabo_p$diabetess==1 & metabo_p$diabetess_doctor==1) | metabo_p$diabetess_medicine==1 | (metabo_p$glucos>= 7 & metabo_p$fasting==0) | (metabo_p$glucos>= 11 & metabo_p$fasting!=0) | (metabo_p$diab_regd <= metabo_p$check_d & !is.na(metabo_p$diab_regd)),1,0)


metabo_sub_chd <- metabo_p[(metabo_p$sub==1  | metabo_p$incchd==1) & metabo_p$fasting!=2,]

## Now read mortality data ##
library(foreign)
dd <- read.dta("/home/andrea/glob/metabo_chd/death_registry_twingene.dta")
dd$deathd <- as.Date(dd$deathd,"%Y-%m-%d")
dd <- dd[,c("twinnr","deathd","dchd")]


metabo_sub2 <- merge(metabo_sub_chd,dd,by="twinnr", all.x=T)


# If missing then assign end of follow-up
metabo_sub2$deathd[is.na(metabo_sub2$deathd)] <- as.Date("20101231","%Y%m%d")
# Survival
metabo_sub2$death_surv <- as.numeric((metabo_sub2$deathd-metabo_sub2$check_d)/365.25)
# If missing death then 0
metabo_sub2$dchd[is.na(metabo_sub2$dchd)] <- 0

# Option 1
# Exclude individuals with non-fatal CHD event #
metabo_sub_nopchdT <- metabo_sub2[metabo_sub_chd$chd==0 | (metabo_sub2$dchd==1 & metabo_sub2$deathd==metabo_sub2$chdd  & !is.na(metabo_sub2$chdd)),]


# Option 2
# Keep individuals with non-fatal CHD event #
metabo_sub_nopchdT <- metabo_sub2


## Log-crp
metabo_sub_nopchdT$logcrp <- log(metabo_sub_nopchdT$crp)
metabo_sub_nopchdT$logtg <- log(metabo_sub_nopchdT$tg)



topfT <- c("M393.231T384.578","M520.340T346.671","M661.528T580.694","M522.356T382.722")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$tg) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi),]


metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$dchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1


## Age in 3 categories ###
metabo_sub_nopchdT4$age3 <- cut(metabo_sub_nopchdT4$age, breaks=c(min(metabo_sub_nopchdT4$age),65,70,max(metabo_sub_nopchdT4$age)),include.lowest = T)

STRATA <- NULL
STRATAM <- NULL
INT <- NULL
MAIN <- NULL
MAIN2 <- NULL
for (i in topfT)
{

	metabo_sub_nopchdT4$x <- scale(as.numeric(metabo_sub_nopchdT4[,i]))
	
	#### CHECK INTERACTION WITH AGE ###
	stratcox2<-coxph(Surv(death_surv,dchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox2,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox2$var)[1],dim(stratcox2$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$death_surv))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox2$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[12])/as.numeric(adjse[12])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(death_surv,dchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

	dfb<-as.matrix(resid(stratcox3,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox3$var)[1],dim(stratcox3$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$death_surv))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox3$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	MAIN <- rbind(MAIN,c(i,paste(round(exp(as.numeric(stratcox3$coefficients[3])),2)," (",round(exp(stratcox3$coefficients[3]+qnorm(0.025)*adjse[3]),2),"-",round(exp(stratcox3$coefficients[3]-qnorm(0.025)*adjse[3]),2),")",sep=""), 2*pnorm(-abs(as.numeric(stratcox3$coefficients[3])/as.numeric(adjse[3])))))
 MAIN2 <- rbind(MAIN2,c(i,as.numeric(stratcox3$coefficients[3]),adjse[3]))


 #### STRATIFIED ####
	stratcox<-coxph(Surv(death_surv,dchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$death_surv))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox$var+gamma


 ###  LOOP OVER DIFFERENT AGE-LEVELS ####
	RES1 <- NULL
	RES2 <- NULL
  for (k in 1:length(levels(metabo_sub_nopchdT4$age3)))
	 {
		newdata0 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=0,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		RES1 <- rbind(RES1,res1)
		res2 <- c(i,lp,se)
		RES2 <- rbind(RES2,res2)
	 }
	
	# Report stratified analysis
	STRATA <- rbind(STRATA,RES1)
	STRATAM <- rbind(STRATAM,RES2[3,])
	
}

## Check if significant interaction
INT
## Check MAIN effect results
MAIN
## Check MAIN effect for meta-analysis
MAIN2
## Check resutls for stratified analysis
STRATA
### Check results for stratified analysis (age > 70), result for meta-analysis
STRATAM


#### Only if interaction use the STRATAM, otherwise just use MAIN ####
#### SO use STRATAM only for M520.340T346.671 ####
REST2 <- rbind(MAIN2[1,],STRATAM[2,],MAIN2[3,],STRATAM[4,])


### META-ANALYSIS ###

library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:length(REST2[,1]))
{
	beta <- c(as.numeric(REST2[i,2]),as.numeric(RESU2[i,2]))
	se <- c(as.numeric(REST2[i,3]),as.numeric(RESU2[i,3]))
	temp <- metagen(beta,se,studlab=c("TwinGene","Ulsam"))
	meta_b <- c(meta_b,temp$TE.random)
	meta_s <- c(meta_s,temp$seTE.random)
}

chd_meta <- cbind(topfT,paste(round(exp(meta_b),2)," (",round(exp(meta_b+qnorm(0.025)*meta_s),2),"-",round(exp(meta_b-qnorm(0.025)*meta_s),2),")",sep=""),2*pnorm(-abs(meta_b/meta_s)))






#################
### NON-FATAL ###
#################


##### ULSAM #####

library(survival)

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/ulsam.metabolomics/data_processed/plasma/ulsam_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]


metabo_sub <- metabo_p[metabo_p$time==0  & metabo_p$double_==0,]

## Now read mortality data ##
library(foreign) 
dd <- read.dta("/home/andrea/glob/metabo_chd/death111013.dta")
dd$p013 <- as.Date(dd$p013,"%Y%m%d")
dd <- dd[,c("pat","chd_death","p013","p002")]

# Merge with mortality data
metabo_sub2 <- merge(metabo_sub,dd,by="pat")

# Indicator if the mortality the incident CHD is due to fatal MI
twexcl <- 	metabo_sub2$pat[metabo_sub2$chd_death==1 & metabo_sub2$p013==metabo_sub2$chdd & !is.na(metabo_sub2$p013) & (metabo_sub2$p013-metabo_sub2$check_d)/365.25 <=10]

metabo_sub_nopchdU <- metabo_sub2[(metabo_sub2$chd==0 | (metabo_sub2$chd==1 & metabo_sub2$chdd > metabo_sub2$check_d)) & !metabo_sub2$pat%in%twexcl,]


#### 10 years follow-up #####
metabo_sub_nopchdU$incchd10 <- ifelse(metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd<=10 & metabo_sub_nopchdU$incchd==1,1,0)
metabo_sub_nopchdU$surv_chd10 <- metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd
metabo_sub_nopchdU$surv_chd10[metabo_sub_nopchdU$surv_chd10>10] <- 10
metabo_sub_nopchdU$surv_chd10[metabo_sub_nopchdU$surv_chd10<0] <- 0.0001
metabo_sub_nopchdU$surv_chd <- metabo_sub_nopchdU$age_exit_chd-metabo_sub_nopchdU$age_entry_chd


topfU <- c("M393.238T397.936","M520.340T357.723","M661.528T590.990","M522.356T395.236")

RESU <- NULL
RESU2 <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU[,i]))+sbp+bmi+smoke+antihyp+ldl+hdl+diab+log(tg), data =metabo_sub_nopchdU)
	res <- c(i,
		paste(round(exp(as.numeric(mod$coefficients[2])),2)," (",
		round(exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),2),"-",
		round(exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]),2),")",sep=""), 	2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])	
	RESU <- rbind(RESU,res)
	RESU2 <- rbind(RESU2,res2)
}





##### TWINGENE #####

library(survival)

load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/twge.metabolomics/data_processed/twge_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]

### Recode diabetes ###
metabo_p$diab_baseline_new <- ifelse((metabo_p$diabetess==1 & metabo_p$diabetess_doctor==1) | metabo_p$diabetess_medicine==1 | (metabo_p$glucos>= 7 & metabo_p$fasting==0) | (metabo_p$glucos>= 11 & metabo_p$fasting!=0) | (metabo_p$diab_regd <= metabo_p$check_d & !is.na(metabo_p$diab_regd)),1,0)


metabo_sub_chd <- metabo_p[(metabo_p$sub==1  | metabo_p$incchd==1) & metabo_p$fasting!=2,]

## Now read mortality data ##
dd <- read.dta("/home/andrea/glob/metabo_chd/death_registry_twingene.dta")
dd$deathd <- as.Date(dd$deathd,"%Y-%m-%d")
dd <- dd[,c("twinnr","deathd","dchd")]


metabo_sub2 <- merge(metabo_sub_chd,dd,by="twinnr", all.x=T)

# Indicator if the mortality the incident CHD is due to fatal MI
twexcl <- 	metabo_sub2$twinnr[metabo_sub2$dchd==1 & metabo_sub2$deathd==metabo_sub2$chdd & !is.na(metabo_sub2$deathd) & !is.na(metabo_sub2$chdd)]

# Exclude fatal MI cases
metabo_sub_nopchdT <- metabo_sub2[(metabo_sub2$chd==0 | metabo_sub2$chd==1 & metabo_sub2$chdd > metabo_sub2$check_d) & !metabo_sub2$twinnr%in%twexcl ,]


## Log-crp
metabo_sub_nopchdT$logcrp <- log(metabo_sub_nopchdT$crp)
metabo_sub_nopchdT$logtg <- log(metabo_sub_nopchdT$tg)



topfT <- c("M393.231T384.578","M520.340T346.671","M661.528T580.694","M522.356T382.722")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$tg) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi),]


metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1


## Age in 3 categories ###
metabo_sub_nopchdT4$age3 <- cut(metabo_sub_nopchdT4$age, breaks=c(min(metabo_sub_nopchdT4$age),65,70,max(metabo_sub_nopchdT4$age)),include.lowest = T)

STRATA <- NULL
STRATAM <- NULL
INT <- NULL
MAIN <- NULL
MAIN2 <- NULL
for (i in topfT)
{

	metabo_sub_nopchdT4$x <- scale(as.numeric(metabo_sub_nopchdT4[,i]))
	
	#### CHECK INTERACTION WITH AGE ###
	stratcox2<-coxph(Surv(surv_chd,incchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox2,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox2$var)[1],dim(stratcox2$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox2$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[12])/as.numeric(adjse[12])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(surv_chd,incchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

	dfb<-as.matrix(resid(stratcox3,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox3$var)[1],dim(stratcox3$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox3$var+gamma
	adjse <- sqrt(diag(adjvar))
	
	MAIN <- rbind(MAIN,c(i,paste(round(exp(as.numeric(stratcox3$coefficients[3])),2)," (",round(exp(stratcox3$coefficients[3]+qnorm(0.025)*adjse[3]),2),"-",round(exp(stratcox3$coefficients[3]-qnorm(0.025)*adjse[3]),2),")",sep=""), 2*pnorm(-abs(as.numeric(stratcox3$coefficients[3])/as.numeric(adjse[3])))))
 MAIN2 <- rbind(MAIN2,c(i,as.numeric(stratcox3$coefficients[3]),adjse[3]))


 #### STRATIFIED ####
	stratcox<-coxph(Surv(surv_chd,incchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
	dfb<-as.matrix(resid(stratcox,type="dfbeta"))

	# Recalculate the correct SE
	strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

	gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
	adjvar<-stratcox$var+gamma


 ###  LOOP OVER DIFFERENT AGE-LEVELS ####
	RES1 <- NULL
	RES2 <- NULL
  for (k in 1:length(levels(metabo_sub_nopchdT4$age3)))
	 {
		newdata0 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=0,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		RES1 <- rbind(RES1,res1)
		res2 <- c(i,lp,se)
		RES2 <- rbind(RES2,res2)
	 }
	
	# Report stratified analysis
	STRATA <- rbind(STRATA,RES1)
	STRATAM <- rbind(STRATAM,RES2[3,])
	
}

## Check if significant interaction
INT
## Check MAIN effect results
MAIN
## Check MAIN effect for meta-analysis
MAIN2
## Check resutls for stratified analysis
STRATA
### Check results for stratified analysis (age > 70), result for meta-analysis
STRATAM


#### Only if interaction use the STRATAM, otherwise just use MAIN ####
#### SO use STRATAM only for M520.340T346.671 ####
REST2 <- rbind(MAIN2[1,],STRATAM[2,],MAIN2[3,],STRATAM[4,])


### META-ANALYSIS ###

library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:length(REST2[,1]))
{
	beta <- c(as.numeric(REST2[i,2]),as.numeric(RESU2[i,2]))
	se <- c(as.numeric(REST2[i,3]),as.numeric(RESU2[i,3]))
	temp <- metagen(beta,se,studlab=c("TwinGene","Ulsam"))
	meta_b <- c(meta_b,temp$TE.random)
	meta_s <- c(meta_s,temp$seTE.random)
}

chd_meta <- cbind(topfT,paste(round(exp(meta_b),2)," (",round(exp(meta_b+qnorm(0.025)*meta_s),2),"-",round(exp(meta_b-qnorm(0.025)*meta_s),2),")",sep=""),2*pnorm(-abs(meta_b/meta_s)))

	

#####################################################################
###### NOT INCLUDED IN THE PAPER - ANALYSIS SPECIFIC ON PIVUS #######
#####################################################################

##### LOAD PIVUS #####
load("/home/andrea/glob/alignment_pivus_small/Results/Final_datasets/step4.Rdata")

## 

# M496.340T365.792 = LysoPC 16:0
# M524.372T425.563 = LysoPC 18:0
# M522.356T384.872 = LysoPC 18:1
# M520.340T347.575 = LysoPC 18:2
# M566.322T352.005 = LysoPC 20:4
# M377.267T388.196 = MG 18:2
# M661.528T579.938 = PE_cer



mPIVUS <- metabo_p[,c("id","M522.356T384.872","M520.340T347.575","M377.267T388.196","M661.528T579.938","check_d","smoke01","diab","tg","age","sex","sbp","diab","tc","hdl","antihyp","ldl","tg")]
cardiphen <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS80_data_CVoutcomes_v2.txt", header=T, sep="\t", stringsAsFactor=F)
mPIVUSF <- merge(mPIVUS,cardiphen,by.x="id", by.y="lpnr")
mPIVUSF$chd <- ifelse(mPIVUSF$mi==1 | mPIVUSF$pci==1, 1, 0)
mPIVUSF$midatestata <- as.Date(mPIVUSF$midatestata, "%d-%b-%y")
mPIVUSF$pcidatestata <- as.Date(mPIVUSF$pcidatestata, "%d-%b-%y")
mPIVUSF$exitdatestata <- as.Date(mPIVUSF$exitdatestata, "%d-%b-%y")
mPIVUSF$chddate <- pmin(mPIVUSF$midatestata,mPIVUSF$pcidatestata)
mPIVUSFi <- mPIVUSF[mPIVUSF$chddate>mPIVUSF$check_d, ]
mPIVUSFi$surv <- mPIVUSFi$chddate-mPIVUSFi$check_d

mod1 <- coxph(Surv(as.numeric(surv),chd)~age+scale(as.numeric(M522.356T384.872))+sex+sbp+diab+ldl+hdl+smoke01+antihyp+log(tg), data= mPIVUSFi)
mod2 <- coxph(Surv(as.numeric(surv),chd)~age+scale(as.numeric(M520.340T347.575))+sex+sbp+diab+ldl+hdl+smoke01+antihyp+log(tg), data= mPIVUSFi)
mod3 <- coxph(Surv(as.numeric(surv),chd)~age+scale(as.numeric(M377.267T388.196))+sex+sbp+diab+ldl+hdl+smoke01+antihyp+log(tg), data= mPIVUSFi)
mod4 <- coxph(Surv(as.numeric(surv),chd)~age+scale(as.numeric(M661.528T579.938))+sex+sbp+diab+ldl+hdl+smoke01+antihyp+log(tg), data= mPIVUSFi)



