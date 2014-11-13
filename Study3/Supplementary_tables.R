#----------------------------------------------
# Filename: Supplementary_tables.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 11JUL2014 - Update the tables names to fit the new paper format, add LDL and TG as covariates, move some tables not in the paper, in the bottom of the program.
# Purpose: Analysis used in the paper. Supplementary table
# Note: 
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small, Pivus small), PIVUS pek tander vs athero6.txt cardiogramplusc4d_data.txt CARDIoGRAM_GWAS_RESULTS.txt SNAPResults.txt
# Data created: 
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/


###############
#### INDEX ####
### 1. SUPPLEMENTARY TABLE 1 - DESCRIPTIVE STATISTICS
### 2. SUPPLEMENTARY TABLE 5 - CLINICAL UTILITY
### 3. SUPPLEMENTARY TABLE 6 - ASSOCIATION BETWEEN MG and incident CHD after adjustment for TG
### 4. SUPPLEMENTARY TABLE 7- ULSAM 20-years FOLLOW-up and top-findings
### 5. SUPPLEMENTARY TABLE 8 - ASSOCIATION WITH CHD ADJUSTING FOR EXTRA COVARIATE
###############
#### LOAD FUNCTIONS ####

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

#C-INDEX
CIND=function(pred,out,time,id,w){
    sumci<- NULL
	sumdi<- NULL
	ci <- NULL
	di <- NULL
for (i in 1:length(id)){
	if (out[i]==1){
		ci <- sum(w[((time[i]>time & pred>pred[i]) | (time>time[i] & pred[i]>pred)) & out==1 & id[i]!=id]) + 2*sum(w[(time>=time[i] & pred[i]>pred) & out==0 & id[i]!=id])
	    di <- sum(w[((time[i]>time & pred<pred[i]) | (time>time[i] & pred[i]<pred)) & out==1 & id[i]!=id]) + 2*sum(w[(time>=time[i] & pred[i]<pred) & out==0 & id[i]!=id])}
	else{
		ci<-0
		di<-0}
	sumci <- sum(c(sumci,ci))
	sumdi <- sum(c(sumdi,di))}
	cindex <- sumci/(sumci+sumdi)
	return(data.frame(cindex))}
	
	

	#NRI
	NRI=function(predwith,predwithout,out,cat,w){
		w <- w[out==0]
		P1 <- predwith
		P2 <- predwithout
		cutP1_e <-cut(P1[out==1], breaks=cat, label=1:(length(cat)-1))
		cutP1_ne <-cut(P1[out==0], breaks=cat, label=1:(length(cat)-1))
		cutP2_e <-cut(P2[out==1], breaks=cat, label=1:(length(cat)-1))
		cutP2_ne <-cut(P2[out==0], breaks=cat, label=1:(length(cat)-1))
		# NON EVENTS
		down_ne <- NULL
		up_ne <- NULL
		dow<-NULL
		up <- NULL
		for (i in 1:length(cutP1_ne)){
			if (as.numeric(cutP2_ne[i]) > as.numeric(cutP1_ne[i])){
				dow = w[i]}
			else {dow=0}
			if  (as.numeric(cutP2_ne[i]) < as.numeric(cutP1_ne[i])){
				up = w[i]}
			else {up=0}
		down_ne <- sum(c(down_ne,dow))
		up_ne <- sum(c(up_ne,up))}
		down_ne <-down_ne/(sum(w))
		up_ne <- up_ne/(sum(w))
		#EVENTS
		down_e <- NULL
		up_e <- NULL
		dow<-NULL
		up <- NULL
		for (i in 1:length(cutP1_e)){
			if (as.numeric(cutP2_e[i]) > as.numeric(cutP1_e[i])){
				dow = 1}
			else {dow=0}
			if  (as.numeric(cutP2_e[i]) < as.numeric(cutP1_e[i])){
				up = 1}
			else {up=0}
		down_e <- sum(c(down_e,dow))
		up_e <- sum(c(up_e,up))}
		down_e <-down_e/sum(out)
		up_e <- up_e/sum(out)
		NRIEXP <- c(((up_e-down_e)+(down_ne-up_ne))*100,(up_e-down_e)*100,(down_ne-up_ne)*100,up_e,down_e,down_ne,up_ne)
		names(NRIEXP) <- c("NRI%","Event NRI%","Non-event NRI%","Up_e","Down_e","down_ne","Up_ne")
	return(NRIEXP)}


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



##### LOAD TWINGENE #####

load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/twge.metabolomics/data_processed/twge_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]

### Recode diabetes ###
metabo_p$diab_baseline_new <- ifelse((metabo_p$diabetess==1 & metabo_p$diabetess_doctor==1) | metabo_p$diabetess_medicine==1 | (metabo_p$glucos>= 7 & metabo_p$fasting==0) | (metabo_p$glucos>= 11 & metabo_p$fasting!=0) | (metabo_p$diab_regd <= metabo_p$check_d & !is.na(metabo_p$diab_regd)),1,0)

metabo_sub_chd <- metabo_p[(metabo_p$sub==1  | metabo_p$incchd==1) & metabo_p$fasting!=2,]
metabo_sub_cvd <- metabo_p[(metabo_p$sub==1  | metabo_p$incchd==1 | metabo_p$incis==1) & metabo_p$fasting!=2,]

metabo_sub_cvd$cvd <- ifelse(metabo_sub_cvd$chd==1 | metabo_sub_cvd$is==1 | metabo_sub_cvd$hs==1,1,0)
metabo_sub_cvd$cvdd <- pmin(metabo_sub_cvd$chdd,metabo_sub_cvd$isd,metabo_sub_cvd$hsd, na.rm=T)
metabo_sub_cvd$inccvd  <- ifelse(metabo_sub_cvd$cvd==1 & metabo_sub_cvd$cvdd > metabo_sub_cvd$check_d,1,0)


# CHD-free subjects
metabo_sub_nopchdT <- metabo_sub_chd[metabo_sub_chd$chd==0 | (metabo_sub_chd$chd==1 & metabo_sub_chd$chdd > metabo_sub_chd$check_d),]

metabo_sub_nopcvdT <- metabo_sub_cvd[metabo_sub_cvd$cvd==0 | (metabo_sub_cvd$cvd==1 & metabo_sub_cvd$cvdd > metabo_sub_cvd$check_d),]
metabo_sub_nopcvdT$surv_cvd <- ifelse(metabo_sub_nopcvdT$inccvd==1, metabo_sub_nopcvdT$cvdd-metabo_sub_nopcvdT$check_d,
	ifelse(metabo_sub_nopcvdT$death==1,metabo_sub_nopcvdT$death_d-metabo_sub_nopcvdT$check_d,as.Date("2010-12-27")-metabo_sub_nopcvdT$check_d))


## Log-crp
metabo_sub_nopchdT$logcrp <- log(metabo_sub_nopchdT$crp)
metabo_sub_nopchdT$logtg <- log(metabo_sub_nopchdT$tg)




##### LOAD PIVUS #####

load("/home/andrea/glob/alignment_pivus_small/Results/Final_datasets/step4.Rdata")
cardiphen <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS pek tander vs athero6.txt", header=T, sep="\t", stringsAsFactor=F)

nnro <- system("awk '{ print NF+1 }' /proj/b2011036/pivus.metabolomics/data_processed/pivus_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]

# Numeric
metabo_p$id <- as.numeric(metabo_p$id)
metabo_chdP <- merge(metabo_p,cardiphen,by="id")
metabo_chdP$logcrp <- log(metabo_chdP$crp.x)
metabo_chdP$logtg <- log(metabo_chdP$tg)





###########################################################
##### 1. SUPPLEMENTARY TABLE 1 - DESCRIPTIVE STATISTICS ###
###########################################################


#### ULSAM ####
RESU <- NULL
for (i in c("age","sbp","bmi","tc","hdl","logcrp","glucose","logtg","ldl"))
{
RESU <- rbind(RESU,c(i,mean(metabo_sub_nopchdU[,i], na.rm=T),sd(metabo_sub_nopchdU[,i],na.rm=T)))
}
median(metabo_sub_nopchdU$surv_chd10)

RESU2 <- NULL
for (i in c("smoke","antihyp","diab","incchd10"))
{
RESU2 <- rbind(RESU2,c(i,table(metabo_sub_nopchdU[,i])[2]/sum(table(metabo_sub_nopchdU[,i]))))
}

quantile(c(metabo_sub_nopchdU$age,metabo_chdP$age), 3/4)
quantile(c(metabo_sub_nopchdU$age,metabo_chdP$age), 1/4)

#### TwinGene ###

library(Hmisc)
library(limma)

metabo_sub_nopchdT$stratum <- ifelse(metabo_sub_nopchdT$incchd==1,0,metabo_sub_nopchdT$stratum)
metabo_sub_nopchdT$weights <- ifelse(metabo_sub_nopchdT$stratum==0,1,
ifelse(metabo_sub_nopchdT$stratum==1,metabo_sub_nopchdT$strata1_n[1]/sum(metabo_sub_nopchdT$stratum==1),
ifelse(metabo_sub_nopchdT$stratum==2,metabo_sub_nopchdT$strata2_n[1]/sum(metabo_sub_nopchdT$stratum==2),
ifelse(metabo_sub_nopchdT$stratum==3,metabo_sub_nopchdT$strata3_n[1]/sum(metabo_sub_nopchdT$stratum==3),metabo_sub_nopchdT$strata4_n[1]/sum(metabo_sub_nopchdT$stratum==4)))))

## Size original cohort ##

metabo_sub_nopchdT$strata1_n[1]+metabo_sub_nopchdT$strata1_n[2]+metabo_sub_nopchdT$strata1_n[3]+metabo_sub_nopchdT$strata1_n[4]

REST <- NULL
for (i in c("age","sbp","bmi","tc","hdl","logcrp","glucos","logtg","ldl"))
{
REST <- rbind(REST,c(i,wtd.mean(metabo_sub_nopchdT[,i], na.rm=T,weights= metabo_sub_nopchdT$weights),sqrt(wtd.var(metabo_sub_nopchdT[,i],na.rm=T,weights= metabo_sub_nopchdT$weights))))
}
weighted.median(metabo_sub_nopchdT$surv_chd, w= metabo_sub_nopchdT$weights)
max(metabo_sub_nopchdT$surv_chd)

REST2 <- NULL
for (i in c("sex","smoke01","antihyp","diab_baseline_new"))
{
REST2 <- rbind(REST2,c(i,wtd.table(metabo_sub_nopchdT[,i],weights=metabo_sub_nopchdT$weights)$ sum.of.weights[2]/sum(wtd.table(metabo_sub_nopchdT[,i],weights=metabo_sub_nopchdT$weights)$ sum.of.weights)))
}

wtd.quantile(c(metabo_sub_nopchdT$age), 3/4, weights=metabo_sub_nopchdT$weights)
wtd.quantile(c(metabo_sub_nopchdT$age), 1/4, weights=metabo_sub_nopchdT$weights)


### PIVUS ###

RESP <- NULL
for (i in c("age","sbp","bmi.x","tc","hdl.x","logcrp","blodsocker","logtg","ldl"))
{
RESP <- rbind(RESP,c(i,mean(metabo_chdP[,i], na.rm=T),sd(metabo_chdP[,i],na.rm=T)))
}

RESP2 <- NULL
for (i in c("sex.x","smoke01","diab","antihyp"))
{
RESP2 <- rbind(RESP2,c(i,table(metabo_chdP[,i])[2]/sum(table(metabo_chdP[,i]))))
}




####################################################
### 2. SUPPLEMENTARY TABLE 5 - CLINICAL UTILITY  ###
####################################################

##### TRADITIONAL APPROACH ####

metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$tc) & !is.na(metabo_sub_nopchdT$smoke) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new),]

metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))


modbasT<-coxph(Surv(surv_chd,incchd) ~ age + sex +sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)


modlysoT<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M520.340T346.671))*age +sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)


modmgT<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M393.231T384.578)) +sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)


modbothT<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M520.340T346.671))*age +  
scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722))*age +
scale(as.numeric(M393.231T384.578)) +
scale(as.numeric(M661.528T580.694)) +
sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)



# Base
cf <- modbasT$coefficient
base_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*metabo_sub_nopchdT4$sbp  + cf[4]*metabo_sub_nopchdT4$smoke + cf[5]*metabo_sub_nopchdT4$antihyp + cf[6]*metabo_sub_nopchdT4$tc + cf[7]*metabo_sub_nopchdT4$hdl +  cf[8]*metabo_sub_nopchdT4$diab_baseline_new

lsurvBT <- survfit(modbasT)
mmBT <- smooth.spline(y=lsurvBT$surv,x=lsurvBT$time/365.25, spar=0.3)
metabo_sub_nopchdT4$base_RF_ORIG <- (1-predict(mmBT,10)$y^exp(base_pred-mean(base_pred)))


# LysoPC
cf <- modlysoT$coefficient
lyso_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671)) + cf[4]*metabo_sub_nopchdT4$sbp  + cf[5]*metabo_sub_nopchdT4$smoke + cf[6]*metabo_sub_nopchdT4$antihyp + cf[7]*metabo_sub_nopchdT4$tc + cf[8]*metabo_sub_nopchdT4$hdl +  cf[9]*metabo_sub_nopchdT4$diab_baseline_new + cf[10]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671))*metabo_sub_nopchdT4$age


lsurvUT <- survfit(modlysoT)
mmUT <- smooth.spline(y=lsurvUT$surv,x=lsurvUT$time/365.25, spar=0.3)

metabo_sub_nopchdT4$lyso_RF_ORIG <- (1-predict(mmUT,10)$y^exp(lyso_pred-mean(lyso_pred)))


# MG 18:2
cf <- modmgT$coefficient
mg_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*scale(as.numeric(metabo_sub_nopchdT4$M393.231T384.578)) + cf[4]*metabo_sub_nopchdT4$sbp  + cf[5]*metabo_sub_nopchdT4$smoke + cf[6]*metabo_sub_nopchdT4$antihyp + cf[7]*metabo_sub_nopchdT4$tc + cf[8]*metabo_sub_nopchdT4$hdl +  cf[9]*metabo_sub_nopchdT4$diab_baseline_new 

lsurvMT <- survfit(modmgT)
mmMT <- smooth.spline(y=lsurvMT$surv,x=lsurvMT$time/365.25, spar=0.3)

metabo_sub_nopchdT4$mg_RF_ORIG <- (1-predict(mmMT,10)$y^exp(mg_pred-mean(mg_pred)))




# BOTH
cf <- modbothT$coefficient
both_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671)) +
cf[4]*scale(as.numeric(metabo_sub_nopchdT4$M522.356T382.722)) +
cf[5]*scale(as.numeric(metabo_sub_nopchdT4$M393.231T384.578)) + 
cf[6]*scale(as.numeric(metabo_sub_nopchdT4$M661.528T580.694)) + 
cf[7]*metabo_sub_nopchdT4$sbp  + cf[8]*metabo_sub_nopchdT4$smoke + cf[9]*metabo_sub_nopchdT4$antihyp + cf[10]*metabo_sub_nopchdT4$tc + cf[11]*metabo_sub_nopchdT4$hdl +  cf[12]*metabo_sub_nopchdT4$diab_baseline_new +
cf[13]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671))*metabo_sub_nopchdT4$age +
cf[14]*scale(as.numeric(metabo_sub_nopchdT4$M522.356T382.722))*metabo_sub_nopchdT4$age

lsurvBOT <- survfit(modbothT)
mmBOT <- smooth.spline(y=lsurvBOT$surv,x=lsurvBOT$time/365.25, spar=0.3)

metabo_sub_nopchdT4$both_RF_ORIG <- (1-predict(mmBOT,10)$y^exp(both_pred-mean(both_pred)))



pdf("test.pdf")
ppBT <- predict(mmBT,c(7,7.5,8,8.5,9,9.5,10))
ppUT <- predict(mmUT,c(7,7.5,8,8.5,9,9.5,10))
ppMT <- predict(mmMT,c(7,7.5,8,8.5,9,9.5,10))
ppBOT <- predict(mmBOT,c(7,7.5,8,8.5,9,9.5,10))

plot(lsurvBT$time/365.25,lsurvBT$surv, ylim=c(0.9,1), xlim=c(0,10))
lines(c(mmBT$x,ppBT$x),c(mmBT$y,ppBT$y),col="red")
lines(c(mmUT$x,ppUT$x),c(mmUT$y,ppUT$y),col="pink")
lines(c(mmMT$x,ppMT$x),c(mmMT$y,ppMT$y),col="grey")
lines(c(mmBOT$x,ppBOT$x),c(mmBOT$y,ppBOT$y),col="yellow")
dev.off()



#### WEIGTHED C-INDEX #####

c_OR_lyso <- CIND(metabo_sub_nopchdT4$lyso_RF_ORIG,metabo_sub_nopchdT4$incchd,metabo_sub_nopchdT4$surv_chd,metabo_sub_nopchdT4$twinnr,metabo_sub_nopchdT4$weights)

c_OR_mg <- CIND(metabo_sub_nopchdT4$mg_RF_ORIG,metabo_sub_nopchdT4$incchd,metabo_sub_nopchdT4$surv_chd,metabo_sub_nopchdT4$twinnr,metabo_sub_nopchdT4$weights)

c_OR_base <-CIND(metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$incchd,metabo_sub_nopchdT4$surv_chd,metabo_sub_nopchdT4$twinnr,metabo_sub_nopchdT4$weights)

c_OR_both <-CIND(metabo_sub_nopchdT4$both_RF_ORIG,metabo_sub_nopchdT4$incchd,metabo_sub_nopchdT4$surv_chd,metabo_sub_nopchdT4$twinnr,metabo_sub_nopchdT4$weights)


### DO A TEST TO COMPARE THE TWO C-INDEXES ###
library(survcomp)
a <- concordance.index(metabo_sub_nopchdT4$both_RF_ORIG,metabo_sub_nopchdT4$surv_chd, metabo_sub_nopchdT4$incchd, weights=metabo_sub_nopchdT4$weights, method="noether")
b <- concordance.index(metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$surv_chd, metabo_sub_nopchdT4$incchd, weights=metabo_sub_nopchdT4$weights, method="noether")
cindex.comp(a, b)


#### WEIGHTED NRI ####

nri_lyso <- NRI(metabo_sub_nopchdT4$lyso_RF_ORIG,metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$incchd,c(0,0.10,0.20,1),metabo_sub_nopchdT4$weights)


nri_mg <- NRI(metabo_sub_nopchdT4$mg_RF_ORIG,metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$incchd,c(0,0.10,0.20,1),metabo_sub_nopchdT4$weights)


nri_both <- NRI(metabo_sub_nopchdT4$both_RF_ORIG,metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$incchd,c(0,0.10,0.20,1),metabo_sub_nopchdT4$weights)




### TwinGene
#### Bootstrapping confidence intervals


NRI_BOTH <- NULL

for (i in 1:200)
{
	set.seed(123+i)
	metabo_sub_nopchdT2 <- metabo_sub_nopchdT[sample(1:nrow(metabo_sub_nopchdT),replace=T),]
	metabo_sub_nopchdT4 <-  metabo_sub_nopchdT2[ !is.na(metabo_sub_nopchdT2$hdl) & !is.na(metabo_sub_nopchdT2$tc) & !is.na(metabo_sub_nopchdT2$smoke) & !is.na(metabo_sub_nopchdT2$sbp) & !is.na(metabo_sub_nopchdT2$antihyp) & !is.na(metabo_sub_nopchdT2$diab_baseline_new),]

	metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

	metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
	ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
	ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
	ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))




	modbasT<-coxph(Surv(surv_chd,incchd) ~ age + sex +sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)

	modbothT<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M520.340T346.671))*age +  
	scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722))*age +
	scale(as.numeric(M393.231T384.578)) +
	scale(as.numeric(M661.528T580.694)) +
	sbp+smoke+antihyp+tc+hdl+diab_baseline_new, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)



	# Base
	cf <- modbasT$coefficient
	base_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*metabo_sub_nopchdT4$sbp  + cf[4]*metabo_sub_nopchdT4$smoke + cf[5]*metabo_sub_nopchdT4$antihyp + cf[6]*metabo_sub_nopchdT4$tc + cf[7]*metabo_sub_nopchdT4$hdl +  cf[8]*metabo_sub_nopchdT4$diab_baseline_new


	lsurvBT <- survfit(modbasT)

	mm <- smooth.spline(y=lsurvBT$surv,x=lsurvBT$time/365.25, spar=0.3)
	pp <- predict(mm,c(7,7.5,8,8.5,9,9.5,10))

	metabo_sub_nopchdT4$base_RF_ORIG <- (1-predict(mm,10)$y^exp(base_pred-mean(base_pred)))


	# BOTH
	cf <- modbothT$coefficient
	both_pred <- cf[1]*metabo_sub_nopchdT4$age + cf[2]*metabo_sub_nopchdT4$sex + cf[3]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671)) +
	cf[4]*scale(as.numeric(metabo_sub_nopchdT4$M522.356T382.722)) +
	cf[5]*scale(as.numeric(metabo_sub_nopchdT4$M393.231T384.578)) + 
	cf[6]*scale(as.numeric(metabo_sub_nopchdT4$M661.528T580.694)) + 
	cf[7]*metabo_sub_nopchdT4$sbp  + cf[8]*metabo_sub_nopchdT4$smoke + cf[9]*metabo_sub_nopchdT4$antihyp + cf[10]*metabo_sub_nopchdT4$tc + cf[11]*metabo_sub_nopchdT4$hdl +  cf[12]*metabo_sub_nopchdT4$diab_baseline_new +
	cf[13]*scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671))*metabo_sub_nopchdT4$age +
	cf[14]*scale(as.numeric(metabo_sub_nopchdT4$M522.356T382.722))*metabo_sub_nopchdT4$age

	lsurvBOT <- survfit(modbothT)
	mmBOT <- smooth.spline(y=lsurvBOT$surv,x=lsurvBOT$time/365.25, spar=0.3)

	metabo_sub_nopchdT4$both_RF_ORIG <- (1-predict(mmBOT,10)$y^exp(both_pred-mean(both_pred)))


	
	nri_both <- NRI(metabo_sub_nopchdT4$both_RF_ORIG,metabo_sub_nopchdT4$base_RF_ORIG,metabo_sub_nopchdT4$incchd,c(0,0.10,0.20,1),metabo_sub_nopchdT4$weights)
	

	NRI_BOTH <- rbind(NRI_BOTH,nri_both)
	print(i)
}

colMeans(NRI_BOTH)
apply(NRI_BOTH,2,function(x) quantile(x,c(0.025)))
apply(NRI_BOTH,2,function(x) quantile(x,c(0.975)))



#### TABLE FOR ESTABLISHED NRI METHOD #####

	P1 <- metabo_sub_nopchdT4$both_RF_ORIG
	P2 <- metabo_sub_nopchdT4$base_RF_ORIG
	out <- metabo_sub_nopchdT4$incchd

	
	cat <- c(0,0.10,0.20,1)
	cutP1_e <-cut(P1[out==1], breaks=cat, label=1:(length(cat)-1))
	cutP1_ne <-cut(P1[out==0], breaks=cat, label=1:(length(cat)-1))
	cutP2_e <-cut(P2[out==1], breaks=cat, label=1:(length(cat)-1))
	cutP2_ne <-cut(P2[out==0], breaks=cat, label=1:(length(cat)-1))
	w <- metabo_sub_nopchdT4$weights
	
	
	wnoout <- w[out==0]
	
	sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==3 & cutP2_ne==3])/(sum(w[out==0]))*100

	
	sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==1 & cutP2_ne==2])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==1 & cutP2_ne==3])/(sum(w[out==0]))*100
	
	
	sum(wnoout[cutP1_ne==2 & cutP2_ne==1])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==2 & cutP2_ne==3])/(sum(w[out==0]))*100
	
	
	
	sum(wnoout[cutP1_ne==3 & cutP2_ne==1])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==3 & cutP2_ne==2])/(sum(w[out==0]))*100
	sum(wnoout[cutP1_ne==3 & cutP2_ne==3])/(sum(w[out==0]))*100
	
	
	
	sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
	sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100
	sum(cutP1_e==3 & cutP2_e==3)/sum(out)*100
	
	
	
	sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
	sum(cutP1_e==1 & cutP2_e==2)/sum(out)*100
	sum(cutP1_e==1 & cutP2_e==3)/sum(out)*100
	
	
	sum(cutP1_e==2 & cutP2_e==1)/sum(out)*100
	sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100
	sum(cutP1_e==2 & cutP2_e==3)/sum(out)*100
	
	
	sum(cutP1_e==3 & cutP2_e==1)/sum(out)*100
	sum(cutP1_e==3 & cutP2_e==2)/sum(out)*100
	sum(cutP1_e==3 & cutP2_e==3)/sum(out)*100
	
	
	

### USING ASCVD
P1 <- metabo_sub_nopchdT4M$both_RF_ASCVD
P2 <- metabo_sub_nopchdT4M$base_RF_ASCVD
out <- metabo_sub_nopchdT4M$incchd


cat <- c(0,0.075,1)
cutP1_e <-cut(P1[out==1], breaks=cat, label=1:(length(cat)-1))
cutP1_ne <-cut(P1[out==0], breaks=cat, label=1:(length(cat)-1))
cutP2_e <-cut(P2[out==1], breaks=cat, label=1:(length(cat)-1))
cutP2_ne <-cut(P2[out==0], breaks=cat, label=1:(length(cat)-1))
w <- metabo_sub_nopchdT4M$weights


wnoout <- w[out==0]

sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==1 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(wnoout[cutP1_ne==2 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100

sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==1 & cutP2_e==2)/sum(out)*100

sum(cutP1_e==2 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100




P1 <- metabo_sub_nopchdT4F$both_RF_ASCVD
P2 <- metabo_sub_nopchdT4F$base_RF_ASCVD
out <- metabo_sub_nopchdT4F$incchd


cat <- c(0,0.075,1)
cutP1_e <-cut(P1[out==1], breaks=cat, label=1:(length(cat)-1))
cutP1_ne <-cut(P1[out==0], breaks=cat, label=1:(length(cat)-1))
cutP2_e <-cut(P2[out==1], breaks=cat, label=1:(length(cat)-1))
cutP2_ne <-cut(P2[out==0], breaks=cat, label=1:(length(cat)-1))
w <- metabo_sub_nopchdT4F$weights


wnoout <- w[out==0]

sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(wnoout[cutP1_ne==1 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==1 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(wnoout[cutP1_ne==2 & cutP2_ne==1])/(sum(w[out==0]))*100
sum(wnoout[cutP1_ne==2 & cutP2_ne==2])/(sum(w[out==0]))*100

sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100

sum(cutP1_e==1 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==1 & cutP2_e==2)/sum(out)*100

sum(cutP1_e==2 & cutP2_e==1)/sum(out)*100
sum(cutP1_e==2 & cutP2_e==2)/sum(out)*100



######################################################################################################
##### 3. SUPPLEMENTARY TABLE 6 - ASSOCIATION BETWEEN MG and incident CHD after adjustment for TG  ####
######################################################################################################


### INCLUDED IN THE SAME MODEL ####

metabo_sub_nopchdT3 <- metabo_sub_nopchdT[!is.na(metabo_sub_nopchdT$logtg),]

metabo_sub_nopchdT3$stratum <- ifelse(metabo_sub_nopchdT3$incchd==1,0,metabo_sub_nopchdT3$stratum)

metabo_sub_nopchdT3$weights <- ifelse(metabo_sub_nopchdT3$stratum==0,1,
ifelse(metabo_sub_nopchdT3$stratum==1,metabo_sub_nopchdT3$strata1_n[1]/sum(metabo_sub_nopchdT3$stratum==1),
ifelse(metabo_sub_nopchdT3$stratum==2,metabo_sub_nopchdT3$strata2_n[1]/sum(metabo_sub_nopchdT3$stratum==2),
ifelse(metabo_sub_nopchdT3$stratum==3,metabo_sub_nopchdT3$strata3_n[1]/sum(metabo_sub_nopchdT3$stratum==3),metabo_sub_nopchdT3$strata4_n[1]/sum(metabo_sub_nopchdT3$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M393.231T384.578))+scale(logtg), data =metabo_sub_nopchdT3,weights=metabo_sub_nopchdT3$weights) 
dfb<-as.matrix(resid(stratcox,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
for (s in 1:4)
	{
			indst<-(1:length(metabo_sub_nopchdT3$surv_chd))[strata_na==s]
			m <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])
			n <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT3$weights)[s+1]))
			if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
	} 
adjvar<-stratcox$var+gamma
adjse <- sqrt(diag(adjvar))
cbind(exp(stratcox$coefficients),adjse,exp(stratcox$coefficients+qnorm(0.025)*adjse),exp(stratcox$coefficients-qnorm(0.025)*adjse), 2*pnorm(-abs(as.numeric(stratcox$coefficients)/as.numeric(adjse))))


#### INCLUDED IN SEPARATE MODELS ####

#### ADJUSTED ####


metabo_sub_nopchdT3 <- metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$logtg),]


metabo_sub_nopchdT3$stratum <- ifelse(metabo_sub_nopchdT3$incchd==1,0,metabo_sub_nopchdT3$stratum)

metabo_sub_nopchdT3$weights <- ifelse(metabo_sub_nopchdT3$stratum==0,1,
ifelse(metabo_sub_nopchdT3$stratum==1,metabo_sub_nopchdT3$strata1_n[1]/sum(metabo_sub_nopchdT3$stratum==1),
ifelse(metabo_sub_nopchdT3$stratum==2,metabo_sub_nopchdT3$strata2_n[1]/sum(metabo_sub_nopchdT3$stratum==2),
ifelse(metabo_sub_nopchdT3$stratum==3,metabo_sub_nopchdT3$strata3_n[1]/sum(metabo_sub_nopchdT3$stratum==3),metabo_sub_nopchdT3$strata4_n[1]/sum(metabo_sub_nopchdT3$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M393.231T384.578))+sbp+smoke+antihyp+ldl+hdl+diab_baseline_new, data =metabo_sub_nopchdT3,weights=metabo_sub_nopchdT3$weights) 
dfb<-as.matrix(resid(stratcox,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
for (s in 1:4)
	{
			indst<-(1:length(metabo_sub_nopchdT3$surv_chd))[strata_na==s]
			m <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])
			n <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT3$weights)[s+1]))
			if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
	} 
adjvar<-stratcox$var+gamma
adjse <- sqrt(diag(adjvar))
cbind(exp(stratcox$coefficients),adjse,exp(stratcox$coefficients+qnorm(0.025)*adjse),exp(stratcox$coefficients-qnorm(0.025)*adjse), 2*pnorm(-abs(as.numeric(stratcox$coefficients)/as.numeric(adjse))))



metabo_sub_nopchdT3 <- metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$logtg),]


metabo_sub_nopchdT3$stratum <- ifelse(metabo_sub_nopchdT3$incchd==1,0,metabo_sub_nopchdT3$stratum)

metabo_sub_nopchdT3$weights <- ifelse(metabo_sub_nopchdT3$stratum==0,1,
ifelse(metabo_sub_nopchdT3$stratum==1,metabo_sub_nopchdT3$strata1_n[1]/sum(metabo_sub_nopchdT3$stratum==1),
ifelse(metabo_sub_nopchdT3$stratum==2,metabo_sub_nopchdT3$strata2_n[1]/sum(metabo_sub_nopchdT3$stratum==2),
ifelse(metabo_sub_nopchdT3$stratum==3,metabo_sub_nopchdT3$strata3_n[1]/sum(metabo_sub_nopchdT3$stratum==3),metabo_sub_nopchdT3$strata4_n[1]/sum(metabo_sub_nopchdT3$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(logtg)+sbp+smoke+antihyp+ldl+hdl+diab_baseline_new, data =metabo_sub_nopchdT3,weights=metabo_sub_nopchdT3$weights) 
dfb<-as.matrix(resid(stratcox,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
for (s in 1:4)
	{
			indst<-(1:length(metabo_sub_nopchdT3$surv_chd))[strata_na==s]
			m <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])
			n <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT3$weights)[s+1]))
			if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
	} 
adjvar<-stratcox$var+gamma
adjse <- sqrt(diag(adjvar))
cbind(exp(stratcox$coefficients),adjse,exp(stratcox$coefficients+qnorm(0.025)*adjse),exp(stratcox$coefficients-qnorm(0.025)*adjse), 2*pnorm(-abs(as.numeric(stratcox$coefficients)/as.numeric(adjse))))




#### UNADJUSTED ####

metabo_sub_nopchdT3 <- metabo_sub_nopchdT[!is.na(metabo_sub_nopchdT$logtg),]


metabo_sub_nopchdT3$stratum <- ifelse(metabo_sub_nopchdT3$incchd==1,0,metabo_sub_nopchdT3$stratum)

metabo_sub_nopchdT3$weights <- ifelse(metabo_sub_nopchdT3$stratum==0,1,
ifelse(metabo_sub_nopchdT3$stratum==1,metabo_sub_nopchdT3$strata1_n[1]/sum(metabo_sub_nopchdT3$stratum==1),
ifelse(metabo_sub_nopchdT3$stratum==2,metabo_sub_nopchdT3$strata2_n[1]/sum(metabo_sub_nopchdT3$stratum==2),
ifelse(metabo_sub_nopchdT3$stratum==3,metabo_sub_nopchdT3$strata3_n[1]/sum(metabo_sub_nopchdT3$stratum==3),metabo_sub_nopchdT3$strata4_n[1]/sum(metabo_sub_nopchdT3$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(M393.231T384.578)), data =metabo_sub_nopchdT3,weights=metabo_sub_nopchdT3$weights) 
dfb<-as.matrix(resid(stratcox,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
for (s in 1:4)
	{
			indst<-(1:length(metabo_sub_nopchdT3$surv_chd))[strata_na==s]
			m <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])
			n <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT3$weights)[s+1]))
			if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
	} 
adjvar<-stratcox$var+gamma
adjse <- sqrt(diag(adjvar))
cbind(exp(stratcox$coefficients),adjse,exp(stratcox$coefficients+qnorm(0.025)*adjse),exp(stratcox$coefficients-qnorm(0.025)*adjse), 2*pnorm(-abs(as.numeric(stratcox$coefficients)/as.numeric(adjse))))



metabo_sub_nopchdT3 <- metabo_sub_nopchdT[!is.na(metabo_sub_nopchdT$logtg),]


metabo_sub_nopchdT3$stratum <- ifelse(metabo_sub_nopchdT3$incchd==1,0,metabo_sub_nopchdT3$stratum)

metabo_sub_nopchdT3$weights <- ifelse(metabo_sub_nopchdT3$stratum==0,1,
ifelse(metabo_sub_nopchdT3$stratum==1,metabo_sub_nopchdT3$strata1_n[1]/sum(metabo_sub_nopchdT3$stratum==1),
ifelse(metabo_sub_nopchdT3$stratum==2,metabo_sub_nopchdT3$strata2_n[1]/sum(metabo_sub_nopchdT3$stratum==2),
ifelse(metabo_sub_nopchdT3$stratum==3,metabo_sub_nopchdT3$strata3_n[1]/sum(metabo_sub_nopchdT3$stratum==3),metabo_sub_nopchdT3$strata4_n[1]/sum(metabo_sub_nopchdT3$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(logtg), data =metabo_sub_nopchdT3,weights=metabo_sub_nopchdT3$weights) 
dfb<-as.matrix(resid(stratcox,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT3$weights))-1

gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
for (s in 1:4)
	{
			indst<-(1:length(metabo_sub_nopchdT3$surv_chd))[strata_na==s]
			m <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])
			n <- as.numeric(table(metabo_sub_nopchdT3$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT3$weights)[s+1]))
			if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
	} 
adjvar<-stratcox$var+gamma
adjse <- sqrt(diag(adjvar))
cbind(exp(stratcox$coefficients),adjse,exp(stratcox$coefficients+qnorm(0.025)*adjse),exp(stratcox$coefficients-qnorm(0.025)*adjse), 2*pnorm(-abs(as.numeric(stratcox$coefficients)/as.numeric(adjse))))





##########################################################################
## 4. SUPPLEMENTARY TABLE  7- ULSAM 20-years FOLLOW-up and top-findings ##
##########################################################################


topfU <- c("M393.238T397.936","M520.340T357.723","M522.356T395.236","M661.528T590.990")
RESU <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(surv_chd20,incchd20) ~ age+scale(as.numeric(metabo_sub_nopchdU[,i]))+sbp+bmi+smoke+antihyp+ldl+hdl+diab+log(tg), data =metabo_sub_nopchdU)
	res <- c(i,exp(as.numeric(mod$coefficients[2])),exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]), 2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))		
	RESU <- rbind(RESU,res)
}

colnames(RESU) <- c("var","HR","HRL","HRH","p-value")





#######################################################################################
##### 5. SUPPLEMENTARY TABLE 8 - ASSOCIATION WITH CHD ADJUSTING FOR EXTRA COVARIATE ###
#######################################################################################



#################################################################
## ADJUST THE ASSOCIATION BETWEEN TOP-FINDINGS AND CHD by CRP ##
################################################################


### ULSAM ###

topfU <- c("M393.238T397.936","M520.340T357.723","M522.356T395.236","M661.528T590.990")
RESUa <- NULL
RESU2 <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU[,i]))+age+sbp+bmi+smoke+antihyp+ldl+hdl+diab+logcrp+log(tg), data =metabo_sub_nopchdU)
	res <- c(i,exp(as.numeric(mod$coefficients[2])),exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]), 2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])		
	RESUa <- rbind(RESUa,res)
	RESU2 <- rbind(RESU2,res2)
}




### TWINGENE ####

topfT <- c("M393.231T384.578","M520.340T346.671","M522.356T382.722","M661.528T580.694")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi) & !is.na(metabo_sub_nopchdT$logcrp) & !is.na(metabo_sub_nopchdT$tg),]


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
	stratcox2<-coxph(Surv(surv_chd,incchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+logcrp+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[13])/as.numeric(adjse[13])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(surv_chd,incchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+logcrp+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

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
	stratcox<-coxph(Surv(surv_chd,incchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+logcrp+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		logcrp=mean(metabo_sub_nopchdT4$logcrp))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		logcrp=mean(metabo_sub_nopchdT4$logcrp))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		res2 <- c(i,lp,se)
		RES1 <- rbind(RES1,res1)
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
REST2 <- rbind(MAIN2[1,],STRATAM[2,],STRATAM[3,],MAIN2[4,])



##### META-ANALYSIS #####

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



####################################################################
## ADJUST THE ASSOCIATION BETWEEN TOP-FINDINGS AND CHD by STATINS ##
####################################################################

### ULSAM ###

## GET STATINS INFORMATION
library (foreign)
d <- read.dta("/home/andrea/glob/alignment_ulsam_small/Data/ingel70+register2002-modified.dta")

d_se <- d[,c("z105","pat")]
 
metabo_sub_nopchdU2 <- merge(metabo_sub_nopchdU,d_se,by="pat")


topfU <- c("M393.238T397.936","M520.340T357.723","M522.356T395.236","M661.528T590.990")
RESUa <- NULL
RESU2 <- NULL
for (i in topfU)
{
	mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU2[,i]))+age+sbp+bmi+smoke+antihyp+ldl+hdl+diab+z105+log(tg), data =metabo_sub_nopchdU2)
	res <- c(i,exp(as.numeric(mod$coefficients[2])),exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]), 2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))	
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])	
	RESUa <- rbind(RESUa,res)
	RESU2 <- rbind(RESU2,res2)
}




#### TWINGENE #####

### Code for statins use ###

metabo_sub_nopchdT4$statins <- ifelse(grepl("Simvastatin",metabo_sub_nopchdT4$hart_medicine_type,ignore.case = T),1,0)


topfT <- c("M393.231T384.578","M520.340T346.671","M522.356T382.722","M661.528T580.694")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi) & !is.na(metabo_sub_nopchdT$tg),]


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
	stratcox2<-coxph(Surv(surv_chd,incchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+statins+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[13])/as.numeric(adjse[13])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(surv_chd,incchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+statins+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

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
	stratcox<-coxph(Surv(surv_chd,incchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+statins+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		statins=mean(metabo_sub_nopchdT4$statins))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		statins=mean(metabo_sub_nopchdT4$statins))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		res2 <- c(i,lp,se)
		RES1 <- rbind(RES1,res1)
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
REST2 <- rbind(MAIN2[1,],STRATAM[2,],STRATAM[3,],MAIN2[4,])

####### META-ANALYSIS #######

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




####################################################################
## ADJUST THE ASSOCIATION BETWEEN TOP-FINDINGS AND CHD by LPPLA2 ##
####################################################################


### TWINGENE ##

### Code LpPLA2###

library (foreign)
d <- read.dta("/home/andrea/glob/alignment_twge_small/Data/lppla2.dta")
d_se <- d[,c("lppla2","twinnr")]

metabo_sub_nopchdTT <- merge(metabo_sub_nopchdT,d_se,by="twinnr")

topfT <- c("M393.231T384.578","M520.340T346.671","M522.356T382.722","M661.528T580.694")


metabo_sub_nopchdT4 <-  metabo_sub_nopchdTT[ !is.na(metabo_sub_nopchdTT$hdl) & !is.na(metabo_sub_nopchdTT$ldl) & !is.na(metabo_sub_nopchdTT$smoke01) & !is.na(metabo_sub_nopchdTT$sbp) & !is.na(metabo_sub_nopchdTT$antihyp) & !is.na(metabo_sub_nopchdTT$diab_baseline_new) & !is.na(metabo_sub_nopchdTT$bmi) & !is.na(metabo_sub_nopchdTT$lppla2) & !is.na(metabo_sub_nopchdTT$tg) ,]


metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

## Age in 3 categories ###
metabo_sub_nopchdT4$age3 <- cut(metabo_sub_nopchdT4$age, breaks=c(min(metabo_sub_nopchdT4$age),65,70,max(metabo_sub_nopchdT4$age)),include.lowest = T)

STRATA <- NULL
INT <- NULL
MAIN <- NULL
for (i in topfT)
{

	metabo_sub_nopchdT4$x <- scale(as.numeric(metabo_sub_nopchdT4[,i]))
	
	#### CHECK INTERACTION WITH AGE ###
	stratcox2<-coxph(Surv(surv_chd,incchd) ~ age + sex +x + x*age +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(lppla2)+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
	
	INT <- rbind(INT,c(i,2*pnorm(-abs(as.numeric(stratcox2$coefficients[13])/as.numeric(adjse[13])))))
	
	

	
	### MAIN EFFECT ####
	stratcox3<-coxph(Surv(surv_chd,incchd) ~ age + sex +x +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(lppla2)+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 

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



 #### STRATIFIED ####
	stratcox<-coxph(Surv(surv_chd,incchd) ~ age3 + sex +x + x*age3 +sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+log(lppla2)+log(tg),  data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights) 
	
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
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		lppla2=mean(log(metabo_sub_nopchdT4$lppla2)))

		newdata1 <- data.frame(
		age3=levels(metabo_sub_nopchdT4$age3)[k],
		x=1,
		sex=mean(metabo_sub_nopchdT4$sex),
		sbp=mean(metabo_sub_nopchdT4$sbp),
		bmi=mean(metabo_sub_nopchdT4$bmi),
		smoke01=mean(metabo_sub_nopchdT4$smoke01),
		antihyp=mean(metabo_sub_nopchdT4$antihyp),
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new),
		lppla2=mean(log(metabo_sub_nopchdT4$lppla2)))

		X0 <- X.coxph(stratcox,newdata=newdata0)
		X1 <- X.coxph(stratcox,newdata=newdata1)
		lp <- predict(stratcox,newdata1) - predict(stratcox,newdata0)
		se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))
		p <- 2*pnorm(-abs(lp/se))

		res1 <- c(i,levels(metabo_sub_nopchdT4$age3)[k],
		paste(round(exp(lp),2)," (",round(exp(lp+qnorm(0.025)*se),2),"-",round(exp(lp-qnorm(0.025)*se),2),")",sep=""),p)
		RES1 <- rbind(RES1,res1)
	 }
	
	# Report stratified analysis
	STRATA <- rbind(STRATA,RES1)
	
}

INT
MAIN
STRATA





#################################################################
### NOT INCLUDED IN THE PAPER - LOOK-UP in  cardiogramplusc4d ###
#################################################################


#### Check if SNPs are in cardiogram ###
a <- read.table("/home/andrea/glob/metabo_chd/cardiogramplusc4d_data.txt", header=T, stringsAsFactor=F)
c <- read.table("/home/andrea/glob/metabo_chd/CARDIoGRAM_GWAS_RESULTS.txt", header=T)

snps <- c("rs75729820","rs8141918","rs174568","rs2048797","rs964184","rs12878001","rs113317091","rs34817779")

a[a[,1]%in%snps,]
c[c[,1]%in%snps,]


## Search for proxy (found with SNAP)
prox <- read.table("/proj/b2011036/nobackup/SNAPResults.txt",sep="\t", header=T, stringsAsFactor=F)

# rs279967 prox for rs75729820, r2=1
# rs1535 prox for rs174568 r2=0.965
# rs10798905 prox for rs113317091, r2=1
# rs7157785 prox for rs12878001, r2=1


snpsp <- c("rs279967","rs1535","rs2048797","rs964184","rs10798905","rs7157785")

## Results for proxy
a[a[,1]%in%snpsp,]
c[c[,1]%in%snpsp,]







###################################################################################
### NOT INCLUDED IN THE PAPER - NRI AND C-STATISTIC IMPROVMENT - ASCVD equation ###
###################################################################################


#### MEN ####

metabo_sub_nopchdT4M <-  metabo_sub_nopcvdT[ !is.na(metabo_sub_nopcvdT$hdl) & !is.na(metabo_sub_nopcvdT$tc) & !is.na(metabo_sub_nopcvdT$smoke) & !is.na(metabo_sub_nopcvdT$sbp) & !is.na(metabo_sub_nopcvdT$antihyp) & !is.na(metabo_sub_nopcvdT$diab_baseline_new) & metabo_sub_nopcvdT$age<80 & metabo_sub_nopcvdT$sex==1,]


metabo_sub_nopchdT4M$stratum <- ifelse(metabo_sub_nopchdT4M$inccvd==1,0,metabo_sub_nopchdT4M$stratum)

metabo_sub_nopchdT4M$weights <- ifelse(metabo_sub_nopchdT4M$stratum==0,1,
ifelse(metabo_sub_nopchdT4M$stratum==1,metabo_sub_nopchdT4M$strata1_n[1]/sum(metabo_sub_nopchdT4M$stratum==1),
ifelse(metabo_sub_nopchdT4M$stratum==2,metabo_sub_nopchdT4M$strata2_n[1]/sum(metabo_sub_nopchdT4M$stratum==2),
ifelse(metabo_sub_nopchdT4M$stratum==3,metabo_sub_nopchdT4M$strata3_n[1]/sum(metabo_sub_nopchdT4M$stratum==3),metabo_sub_nopchdT4M$strata4_n[1]/sum(metabo_sub_nopchdT4M$stratum==4)))))


base_LPM <- 12.344*log(metabo_sub_nopchdT4M$age) + 
11.853*log(metabo_sub_nopchdT4M$tc*38.66976)  + 
-7.990*log(metabo_sub_nopchdT4M$hdl*38.66976) + 
1.764*log(metabo_sub_nopchdT4M$sbp) + 
7.837*metabo_sub_nopchdT4M$smoke + 
0.658*metabo_sub_nopchdT4M$diab_baseline_new +
-2.664*(log(metabo_sub_nopchdT4M$tc*38.66976)*log(metabo_sub_nopchdT4M$age)) + 
1.769*(log(metabo_sub_nopchdT4M$hdl*38.66976)*log(metabo_sub_nopchdT4M$age)) + 
-1.795*(metabo_sub_nopchdT4M$smoke*log(metabo_sub_nopchdT4M$age))


## Estimate the coefficients
modASCVDM<-coxph(Surv(surv_cvd,inccvd) ~ scale(as.numeric(M520.340T346.671)) + offset(scale(base_LPM)) + scale(as.numeric(M520.340T346.671))*age, data =metabo_sub_nopchdT4M,weights=metabo_sub_nopchdT4M$weights)

modASCVDM2<-coxph(Surv(surv_cvd,inccvd) ~ scale(as.numeric(M393.231T384.578)) + offset(scale(base_LPM)), data =metabo_sub_nopchdT4M,weights=metabo_sub_nopchdT4M$weights)

modASCVDM3<-coxph(Surv(surv_cvd,inccvd) ~  scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M393.231T384.578)) + scale(as.numeric(M520.340T346.671)):age + scale(as.numeric(M661.528T580.694)) + scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722)):age + offset(scale(base_LPM)), data =metabo_sub_nopchdT4M,weights=metabo_sub_nopchdT4M$weights)






# Add new coefficients to the existing risk estimator
lyso_LPM <-base_LPM + modASCVDM$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))+modASCVDM$coefficients[3]*(scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))*metabo_sub_nopchdT4M$age)

mg_LPM <-base_LPM + modASCVDM2$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4M$M393.231T384.578))

both_LPM <-base_LPM + 
modASCVDM3$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))+
modASCVDM3$coefficients[2]*scale(as.numeric(metabo_sub_nopchdT4M$M393.231T384.578)) +
modASCVDM3$coefficients[3]*scale(as.numeric(metabo_sub_nopchdT4M$M661.528T580.694)) +
modASCVDM3$coefficients[4]*scale(as.numeric(metabo_sub_nopchdT4M$M522.356T382.722)) +
modASCVDM3$coefficients[5]*(scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))*metabo_sub_nopchdT4M$age) +
modASCVDM3$coefficients[6]*(scale(as.numeric(metabo_sub_nopchdT4M$M522.356T382.722))*metabo_sub_nopchdT4M$age)


lsurvM <- survfit(Surv(surv_cvd/365.25, inccvd) ~ 1, data=metabo_sub_nopchdT4M, type='fleming',weights=metabo_sub_nopchdT4M$weights)


mm <- lm(lsurvM$surv~lsurvM$time)
pp <- mm$coefficients[1]+c(7,7.5,8,8.5,9,9.5,10)*mm$coefficients[2]

pdf("test.pdf")
plot(lsurvM$time,lsurvM$surv, ylim=c(0.9,1), xlim=c(0,10))
lines(c(lsurvM$time,c(7,7.5,8,8.5,9,9.5,10)),c(lsurvM$surv,pp),col="red")
dev.off()

metabo_sub_nopchdT4M$base_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(base_LPM-mean(base_LPM)))

metabo_sub_nopchdT4M$lyso_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(lyso_LPM-mean(lyso_LPM)))

metabo_sub_nopchdT4M$mg_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(mg_LPM-mean(mg_LPM)))

metabo_sub_nopchdT4M$both_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(both_LPM-mean(both_LPM)))



#### WEIGHTED C-INDEX ####
c_ASCVD_lysoM <- CIND(metabo_sub_nopchdT4M$lyso_RF_ASCVD,metabo_sub_nopchdT4M$incchd,metabo_sub_nopchdT4M$surv_cvd,metabo_sub_nopchdT4M$twinnr,metabo_sub_nopchdT4M$weights)
c_ASCVD_baseM <-CIND(metabo_sub_nopchdT4M$base_RF_ASCVD,metabo_sub_nopchdT4M$incchd,metabo_sub_nopchdT4M$surv_cvd,metabo_sub_nopchdT4M$twinnr,metabo_sub_nopchdT4M$weights)

c_ASCVD_mgM <- CIND(metabo_sub_nopchdT4M$mg_RF_ASCVD,metabo_sub_nopchdT4M$incchd,metabo_sub_nopchdT4M$surv_cvd,metabo_sub_nopchdT4M$twinnr,metabo_sub_nopchdT4M$weights)


c_ASCVD_bothM <- CIND(metabo_sub_nopchdT4M$both_RF_ASCVD,metabo_sub_nopchdT4M$incchd,metabo_sub_nopchdT4M$surv_cvd,metabo_sub_nopchdT4M$twinnr,metabo_sub_nopchdT4M$weights)




#### WEIGHTED NRI ####
nri_lysoM <- NRI(metabo_sub_nopchdT4M$lyso_RF_ASCVD,metabo_sub_nopchdT4M$base_RF_ASCVD,metabo_sub_nopchdT4M$inccvd,c(0,0.075,1),metabo_sub_nopchdT4M$weights)


nri_mgM <- NRI(metabo_sub_nopchdT4M$mg_RF_ASCVD,metabo_sub_nopchdT4M$base_RF_ASCVD,metabo_sub_nopchdT4M$inccvd,c(0,0.075,1),metabo_sub_nopchdT4M$weights)


nri_bothM <- NRI(metabo_sub_nopchdT4M$both_RF_ASCVD,metabo_sub_nopchdT4M$base_RF_ASCVD,metabo_sub_nopchdT4M$inccvd,c(0,0.075,1),metabo_sub_nopchdT4M$weights)





#### WOMEN ####



metabo_sub_nopchdT4F <-  metabo_sub_nopcvdT[ !is.na(metabo_sub_nopcvdT$hdl) & !is.na(metabo_sub_nopcvdT$tc) & !is.na(metabo_sub_nopcvdT$smoke) & !is.na(metabo_sub_nopcvdT$sbp) & !is.na(metabo_sub_nopcvdT$antihyp) & !is.na(metabo_sub_nopcvdT$diab_baseline_new) & metabo_sub_nopcvdT$age<80 & metabo_sub_nopcvdT$sex==2,]


metabo_sub_nopchdT4F$stratum <- ifelse(metabo_sub_nopchdT4F$inccvd==1,0,metabo_sub_nopchdT4F$stratum)

metabo_sub_nopchdT4F$weights <- ifelse(metabo_sub_nopchdT4F$stratum==0,1,
ifelse(metabo_sub_nopchdT4F$stratum==1,metabo_sub_nopchdT4F$strata1_n[1]/sum(metabo_sub_nopchdT4F$stratum==1),
ifelse(metabo_sub_nopchdT4F$stratum==2,metabo_sub_nopchdT4F$strata2_n[1]/sum(metabo_sub_nopchdT4F$stratum==2),
ifelse(metabo_sub_nopchdT4F$stratum==3,metabo_sub_nopchdT4F$strata3_n[1]/sum(metabo_sub_nopchdT4F$stratum==3),metabo_sub_nopchdT4F$strata4_n[1]/sum(metabo_sub_nopchdT4F$stratum==4)))))






base_LPF <- -29.799*log(metabo_sub_nopchdT4F$age) + 
4.884*I(log(metabo_sub_nopchdT4F$age)^2) + 
13.540*log(metabo_sub_nopchdT4F$tc*38.66976)  + 
-13.578*log(metabo_sub_nopchdT4F$hdl*38.66976) + 
2.091*log(metabo_sub_nopchdT4F$sbp) + 
7.574*metabo_sub_nopchdT4F$smoke + 
0.661*metabo_sub_nopchdT4F$diab_baseline_new +
-3.114*(log(metabo_sub_nopchdT4F$tc*38.66976)*log(metabo_sub_nopchdT4F$age)) + 
3.149*(log(metabo_sub_nopchdT4F$hdl*38.66976)*log(metabo_sub_nopchdT4F$age)) + 
-1.665*(metabo_sub_nopchdT4F$smoke*log(metabo_sub_nopchdT4F$age))



## Estimate the coefficients
modASCVDF<-coxph(Surv(surv_cvd,inccvd) ~ scale(as.numeric(M520.340T346.671)) + offset(scale(base_LPF)) + scale(as.numeric(M520.340T346.671))*age, data =metabo_sub_nopchdT4F,weights=metabo_sub_nopchdT4F$weights)

modASCVDF2<-coxph(Surv(surv_cvd,inccvd) ~ scale(as.numeric(M393.231T384.578)) + offset(scale(base_LPF)), data =metabo_sub_nopchdT4F,weights=metabo_sub_nopchdT4F$weights)

modASCVDF3<-coxph(Surv(surv_cvd,inccvd) ~  scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M393.231T384.578)) + scale(as.numeric(M520.340T346.671)):age + scale(as.numeric(M661.528T580.694)) + scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722)):age + offset(scale(base_LPF)), data =metabo_sub_nopchdT4F,weights=metabo_sub_nopchdT4F$weights)



lyso_LPF <-base_LPF + modASCVDF$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671)) + modASCVDF$coefficients[3]*(scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671))*metabo_sub_nopchdT4F$age)


mg_LPF <-base_LPF + modASCVDF2$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4F$M393.231T384.578))

both_LPF <-base_LPF + 
modASCVDF3$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671))+
modASCVDF3$coefficients[2]*scale(as.numeric(metabo_sub_nopchdT4F$M393.231T384.578)) +
modASCVDF3$coefficients[3]*scale(as.numeric(metabo_sub_nopchdT4F$M661.528T580.694)) +
modASCVDF3$coefficients[4]*scale(as.numeric(metabo_sub_nopchdT4F$M522.356T382.722)) +
modASCVDF3$coefficients[5]*(scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671))*metabo_sub_nopchdT4F$age) +
modASCVDF3$coefficients[6]*(scale(as.numeric(metabo_sub_nopchdT4F$M522.356T382.722))*metabo_sub_nopchdT4F$age)




lsurvF <- survfit(Surv(surv_cvd/365.25, inccvd) ~ 1, data=metabo_sub_nopchdT4F, type='fleming',weights=metabo_sub_nopchdT4F$weights)


mm <- lm(lsurvF$surv~lsurvF$time)
pp <- mm$coefficients[1]+c(7,7.5,8,8.5,9,9.5,10)*mm$coefficients[2]

pdf("test.pdf")
plot(lsurvF$time,lsurvF$surv, ylim=c(0.9,1), xlim=c(0,10))
lines(c(lsurvF$time,c(7,7.5,8,8.5,9,9.5,10)),c(lsurvF$surv,pp),col="red")
dev.off()

metabo_sub_nopchdT4F$base_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(base_LPF-mean(base_LPF)))

metabo_sub_nopchdT4F$lyso_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(lyso_LPF-mean(lyso_LPF)))

metabo_sub_nopchdT4F$mg_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(mg_LPF-mean(mg_LPF)))

metabo_sub_nopchdT4F$both_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(both_LPF-mean(both_LPF)))



#### WEIGHTED C-INDEX ####
c_ASCVD_lysoF <- CIND(metabo_sub_nopchdT4F$lyso_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,metabo_sub_nopchdT4F$surv_cvd,metabo_sub_nopchdT4F$twinnr,metabo_sub_nopchdT4F$weights)
c_ASCVD_baseF <-CIND(metabo_sub_nopchdT4F$base_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,metabo_sub_nopchdT4F$surv_cvd,metabo_sub_nopchdT4F$twinnr,metabo_sub_nopchdT4F$weights)

c_ASCVD_mgF <- CIND(metabo_sub_nopchdT4F$mg_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,metabo_sub_nopchdT4F$surv_cvd,metabo_sub_nopchdT4F$twinnr,metabo_sub_nopchdT4F$weights)


c_ASCVD_bothF <- CIND(metabo_sub_nopchdT4F$both_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,metabo_sub_nopchdT4F$surv_cvd,metabo_sub_nopchdT4F$twinnr,metabo_sub_nopchdT4F$weights)




#### WEIGHTED NRI ####
nri_lysoF <- NRI(metabo_sub_nopchdT4F$lyso_RF_ASCVD,metabo_sub_nopchdT4F$base_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,c(0,0.075,1),metabo_sub_nopchdT4F$weights)


nri_mgF <- NRI(metabo_sub_nopchdT4F$mg_RF_ASCVD,metabo_sub_nopchdT4F$base_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,c(0,0.075,1),metabo_sub_nopchdT4F$weights)

nri_bothF <- NRI(metabo_sub_nopchdT4F$both_RF_ASCVD,metabo_sub_nopchdT4F$base_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,c(0,0.075,1),metabo_sub_nopchdT4F$weights)



#### BOOTSTRAP #####


NRI_BOTHM <- NULL
NRI_BOTHF <- NULL

for (i in 1:200)
{
	set.seed(123+i)
	metabo_sub_nopcvdT2 <- metabo_sub_nopcvdT[sample(1:nrow(metabo_sub_nopcvdT),replace=T),]

	metabo_sub_nopchdT4M <-  metabo_sub_nopcvdT2[ !is.na(metabo_sub_nopcvdT2$hdl) & !is.na(metabo_sub_nopcvdT2$tc) & !is.na(metabo_sub_nopcvdT2$smoke) & !is.na(metabo_sub_nopcvdT2$sbp) & !is.na(metabo_sub_nopcvdT2$antihyp) & !is.na(metabo_sub_nopcvdT2$diab_baseline_new) & metabo_sub_nopcvdT2$age<80 & metabo_sub_nopcvdT2$sex==1,]


	metabo_sub_nopchdT4M$stratum <- ifelse(metabo_sub_nopchdT4M$inccvd==1,0,metabo_sub_nopchdT4M$stratum)

	metabo_sub_nopchdT4M$weights <- ifelse(metabo_sub_nopchdT4M$stratum==0,1,
	ifelse(metabo_sub_nopchdT4M$stratum==1,metabo_sub_nopchdT4M$strata1_n[1]/sum(metabo_sub_nopchdT4M$stratum==1),
	ifelse(metabo_sub_nopchdT4M$stratum==2,metabo_sub_nopchdT4M$strata2_n[1]/sum(metabo_sub_nopchdT4M$stratum==2),
	ifelse(metabo_sub_nopchdT4M$stratum==3,metabo_sub_nopchdT4M$strata3_n[1]/sum(metabo_sub_nopchdT4M$stratum==3),metabo_sub_nopchdT4M$strata4_n[1]/sum(metabo_sub_nopchdT4M$stratum==4)))))


	base_LPM <- 12.344*log(metabo_sub_nopchdT4M$age) + 
	11.853*log(metabo_sub_nopchdT4M$tc*38.66976)  + 
	-7.990*log(metabo_sub_nopchdT4M$hdl*38.66976) + 
	1.764*log(metabo_sub_nopchdT4M$sbp) + 
	7.837*metabo_sub_nopchdT4M$smoke + 
	0.658*metabo_sub_nopchdT4M$diab_baseline_new +
	-2.664*(log(metabo_sub_nopchdT4M$tc*38.66976)*log(metabo_sub_nopchdT4M$age)) + 
	1.769*(log(metabo_sub_nopchdT4M$hdl*38.66976)*log(metabo_sub_nopchdT4M$age)) + 
	-1.795*(metabo_sub_nopchdT4M$smoke*log(metabo_sub_nopchdT4M$age))


	modASCVDM3<-coxph(Surv(surv_cvd,inccvd) ~  scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M393.231T384.578)) + scale(as.numeric(M520.340T346.671)):age + scale(as.numeric(M661.528T580.694)) + scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722)):age + offset(scale(base_LPM)), data =metabo_sub_nopchdT4M,weights=metabo_sub_nopchdT4M$weights)

	both_LPM <-base_LPM + 
	modASCVDM3$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))+
	modASCVDM3$coefficients[2]*scale(as.numeric(metabo_sub_nopchdT4M$M393.231T384.578)) +
	modASCVDM3$coefficients[3]*scale(as.numeric(metabo_sub_nopchdT4M$M661.528T580.694)) +
	modASCVDM3$coefficients[4]*scale(as.numeric(metabo_sub_nopchdT4M$M522.356T382.722)) +
	modASCVDM3$coefficients[5]*(scale(as.numeric(metabo_sub_nopchdT4M$M520.340T346.671))*metabo_sub_nopchdT4M$age) +
	modASCVDM3$coefficients[6]*(scale(as.numeric(metabo_sub_nopchdT4M$M522.356T382.722))*metabo_sub_nopchdT4M$age)


	lsurvM <- survfit(Surv(surv_cvd/365.25, inccvd) ~ 1, data=metabo_sub_nopchdT4M, type='fleming',weights=metabo_sub_nopchdT4M$weights)

	mm <- lm(lsurvM$surv~lsurvM$time)
	pp <- mm$coefficients[1]+c(7,7.5,8,8.5,9,9.5,10)*mm$coefficients[2]

	metabo_sub_nopchdT4M$base_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(base_LPM-mean(base_LPM)))

	metabo_sub_nopchdT4M$both_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(both_LPM-mean(both_LPM)))


	nri_bothM <- NRI(metabo_sub_nopchdT4M$both_RF_ASCVD,metabo_sub_nopchdT4M$base_RF_ASCVD,metabo_sub_nopchdT4M$inccvd,c(0,0.075,1),metabo_sub_nopchdT4M$weights)

	
	NRI_BOTHM <- rbind(NRI_BOTHM,nri_bothM)
	
	
	
	metabo_sub_nopchdT4F <-  metabo_sub_nopcvdT2[ !is.na(metabo_sub_nopcvdT2$hdl) & !is.na(metabo_sub_nopcvdT2$tc) & !is.na(metabo_sub_nopcvdT2$smoke) & !is.na(metabo_sub_nopcvdT2$sbp) & !is.na(metabo_sub_nopcvdT2$antihyp) & !is.na(metabo_sub_nopcvdT2$diab_baseline_new) & metabo_sub_nopcvdT2$age<80 & metabo_sub_nopcvdT2$sex==2,]


	metabo_sub_nopchdT4F$stratum <- ifelse(metabo_sub_nopchdT4F$inccvd==1,0,metabo_sub_nopchdT4F$stratum)

	metabo_sub_nopchdT4F$weights <- ifelse(metabo_sub_nopchdT4F$stratum==0,1,
	ifelse(metabo_sub_nopchdT4F$stratum==1,metabo_sub_nopchdT4F$strata1_n[1]/sum(metabo_sub_nopchdT4F$stratum==1),
	ifelse(metabo_sub_nopchdT4F$stratum==2,metabo_sub_nopchdT4F$strata2_n[1]/sum(metabo_sub_nopchdT4F$stratum==2),
	ifelse(metabo_sub_nopchdT4F$stratum==3,metabo_sub_nopchdT4F$strata3_n[1]/sum(metabo_sub_nopchdT4F$stratum==3),metabo_sub_nopchdT4F$strata4_n[1]/sum(metabo_sub_nopchdT4F$stratum==4)))))



	base_LPF <- -29.799*log(metabo_sub_nopchdT4F$age) + 
	4.884*I(log(metabo_sub_nopchdT4F$age)^2) + 
	13.540*log(metabo_sub_nopchdT4F$tc*38.66976)  + 
	-13.578*log(metabo_sub_nopchdT4F$hdl*38.66976) + 
	2.091*log(metabo_sub_nopchdT4F$sbp) + 
	7.574*metabo_sub_nopchdT4F$smoke + 
	0.661*metabo_sub_nopchdT4F$diab_baseline_new +
	-3.114*(log(metabo_sub_nopchdT4F$tc*38.66976)*log(metabo_sub_nopchdT4F$age)) + 
	3.149*(log(metabo_sub_nopchdT4F$hdl*38.66976)*log(metabo_sub_nopchdT4F$age)) + 
	-1.665*(metabo_sub_nopchdT4F$smoke*log(metabo_sub_nopchdT4F$age))


	modASCVDF3<-coxph(Surv(surv_cvd,inccvd) ~  scale(as.numeric(M520.340T346.671)) + scale(as.numeric(M393.231T384.578)) + scale(as.numeric(M520.340T346.671)):age + scale(as.numeric(M661.528T580.694)) + scale(as.numeric(M522.356T382.722)) + scale(as.numeric(M522.356T382.722)):age + offset(scale(base_LPF)), data =metabo_sub_nopchdT4F,weights=metabo_sub_nopchdT4F$weights)


	both_LPF <-base_LPF + 
	modASCVDF3$coefficients[1]*scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671))+
	modASCVDF3$coefficients[2]*scale(as.numeric(metabo_sub_nopchdT4F$M393.231T384.578)) +
	modASCVDF3$coefficients[3]*scale(as.numeric(metabo_sub_nopchdT4F$M661.528T580.694)) +
	modASCVDF3$coefficients[4]*scale(as.numeric(metabo_sub_nopchdT4F$M522.356T382.722)) +
	modASCVDF3$coefficients[5]*(scale(as.numeric(metabo_sub_nopchdT4F$M520.340T346.671))*metabo_sub_nopchdT4F$age) +
	modASCVDF3$coefficients[6]*(scale(as.numeric(metabo_sub_nopchdT4F$M522.356T382.722))*metabo_sub_nopchdT4F$age)


	lsurvF <- survfit(Surv(surv_cvd/365.25, inccvd) ~ 1, data=metabo_sub_nopchdT4F, type='fleming',weights=metabo_sub_nopchdT4F$weights)


	mm <- lm(lsurvF$surv~lsurvF$time)
	pp <- mm$coefficients[1]+c(7,7.5,8,8.5,9,9.5,10)*mm$coefficients[2]


	metabo_sub_nopchdT4F$base_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(base_LPF-mean(base_LPF)))

	metabo_sub_nopchdT4F$both_RF_ASCVD <- (1-(mm$coefficients[1]+10*mm$coefficients[2])^exp(both_LPF-mean(both_LPF)))

	nri_bothF <- NRI(metabo_sub_nopchdT4F$both_RF_ASCVD,metabo_sub_nopchdT4F$base_RF_ASCVD,metabo_sub_nopchdT4F$inccvd,c(0,0.075,1),metabo_sub_nopchdT4F$weights)
	
	NRI_BOTHF <- rbind(NRI_BOTHF,nri_bothF)
		
	
	print(i)
}


colMeans(NRI_BOTHF)
apply(NRI_BOTHF,2,function(x) quantile(x,c(0.025)))
apply(NRI_BOTHF,2,function(x) quantile(x,c(0.975)))

colMeans(NRI_BOTHM)
apply(NRI_BOTHM,2,function(x) quantile(x,c(0.025)))
apply(NRI_BOTHM,2,function(x) quantile(x,c(0.975)))




