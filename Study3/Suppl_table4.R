#----------------------------------------------
# Filename: Suppl_table4.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 11JUL2014 - Rename the file to reflect new paper format and include LDL cholesterol and TG as covariates
# Purpose: Analysis used in the paper. Supplementary table 4.
# Note: Association between LysoPCs and LysoPCs ratios and incident CHD
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small, Pivus small)
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


## TWINGENE ##

load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## 
# M494.325T325.856 = LysoPC 16:1
# M496.340T364.064 = LysoPC 16:0
# M522.356T382.722 = LysoPC 18:1
# M524.372T423.532 = LysoPC 18:0
# M520.340T346.671 = LysoPC 18:2
# M566.322T340.659 = LysoPC 20:4

# M730.538T619.424 = PC(32:2) 
# M732.554T640.164 = PC(32:1)
# M734.560T640.164 = PC(32:0)
# M782.570T652.163 = PC(36:4)
# M812.616T687.249 = PC(38:3)


# M760.585T667.863 = PC(34:1)
# M758.570T649.594 = PC(34:2)
# M788.617T694.035 = PC(36:1)
# M786.601T677.527 = PC(36:2)


# M285.279T515.220 = STEARIC ACID
# M283.264T452.573 = OLEIC / ELAIDIC ACID
# M321.270T463.343 = VACCENIC ACID
# M551.434T435.908 = PALMITIC ACID
# M293.179T392.421 = PALMITOLEIC
# M281.010T31.737 = GLYCEROL PC
# M319.195T415.326 = LIONELIC ACID 

# M118.087T33.403 = Betaine
# M105.111T32.168 = Choline


mTWINGENEF <- metabo_p[,c("twinnr","M494.325T325.856","M496.340T364.064","M522.356T382.722","M524.372T423.532","M520.340T346.671","M566.322T340.659",
"M730.538T619.424","M732.554T640.164","M734.560T640.164","M782.570T652.163",
"M786.601T677.527","M812.616T687.249",
"M285.279T515.220","M283.264T452.573","M321.270T463.343","M551.434T435.908","M293.179T392.421","M281.010T31.737","M319.195T415.326",
"M118.087T33.403","M105.111T32.168",
"M760.585T667.863","M758.570T649.594","M788.617T694.035",colnames(metabo_p)[9756:ncol(metabo_p)])]

colnames(mTWINGENEF)[2:25] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC32_2","PC32_1","PC32_0","PC36_4","PC36_2","PC38_3",
"STEARIC","OLEIC","VACCENIC","PALMITIC","PALMITOLEIC","GLYPC","LIONELIC",
"BETAINE","CHOLINE",
"PC34_1","PC34_2","PC36_1")

mTWINGENEF$logcrp <- log(mTWINGENEF$crp)
mTWINGENEF$logtg <- log(mTWINGENEF$tg)



## ULSAM ##

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## 
# M494.324T324.997 = LysoPC 16:1
# M496.340T375.465 = LysoPC 16:0
# M522.356T395.236 = LysoPC 18:1
# M524.372T436.589 = LysoPC 18:0
# M520.340T357.723 = LysoPC 18:2
# M566.322T351.269 = LysoPC 20:4

# M730.538T628.356 = PC(32:2) 
# M732.554T648.330 = PC(32:1)
# M734.560T648.335 = PC(32:0)
# M782.575T647.537 = PC(36:4)
# M812.616T693.845 = PC(38:3) 


# M760.594T674.770 = PC(34:1)
# M758.569T657.575 = PC(34:2)
# M788.616T700.812 = PC(36:1)
# M786.609T683.883 = PC(36:2)


# M285.280T528.734 = STEARIC ACID
# M283.264T466.637 = OLEIC / ELAIDIC ACID
# M321.270T477.975 = VACCENIN ACID
# M551.434T449.611 = PALMITIC ACID
# M293.179T405.566 = PALMITOLEIC ACID
# M281.010T31.036 = GLYCEROL PC
# M319.195T428.855 = LIONELIC ACID

# M118.087T34.235 = Betaine
# M105.111T32.235 = Choline




mULSAMF <- metabo_p[metabo_p$time==0 & metabo_p$double_==0,c("pat_time","M494.324T324.997","M496.340T375.465","M522.356T395.236","M524.372T436.589","M520.340T357.723","M566.322T351.269",
"M730.538T628.356","M732.554T648.330","M734.560T648.335","M782.575T647.537","M786.609T683.883","M812.616T693.845",
"M285.280T528.734","M283.264T466.637","M321.270T477.975","M551.434T449.611","M293.179T405.566","M281.010T31.036","M319.195T428.855",
"M118.087T34.235","M105.111T32.235",
"M760.594T674.770","M758.569T657.575","M788.616T700.812",colnames(metabo_p)[10163:ncol(metabo_p)])]

colnames(mULSAMF)[2:25] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC32_2","PC32_1","PC32_0","PC36_4","PC36_2","PC38_3",
"STEARIC","OLEIC","VACCENIC","PALMITIC","PALMITOLEIC","GLYPC","LIONELIC",
"BETAINE","CHOLINE",
"PC34_1","PC34_2","PC36_1")

mULSAMF$logcrp <- log(mULSAMF$crp)
mULSAMF$logtg <- log(mULSAMF$tg)



###############
#### ULSAM ####
###############


##### SELECT CHD EVENTS #####
mULSAMFC <- mULSAMF[mULSAMF$chd==0 | (mULSAMF$chd==1 & mULSAMF$chdd > mULSAMF$check_d),]

#### 10 years follow-up #####
mULSAMFC$incchd10 <- ifelse(mULSAMFC$age_exit_chd-mULSAMFC$age_entry_chd<=10 & mULSAMFC$incchd==1,1,0)
mULSAMFC$surv_chd10 <- mULSAMFC$age_exit_chd-mULSAMFC$age_entry_chd
mULSAMFC$surv_chd10[mULSAMFC$surv_chd10>10] <- 10
mULSAMFC$surv_chd10[mULSAMFC$surv_chd10<0] <- 0.0001
mULSAMFC$surv_chd <- mULSAMFC$age_exit_chd-mULSAMFC$age_entry_chd


#### 20 years follow-up #####
mULSAMFC$incchd20 <- ifelse(mULSAMFC$age_exit_chd-mULSAMFC$age_entry_chd<=20 & mULSAMFC$incchd==1,1,0)
mULSAMFC$surv_chd20 <- mULSAMFC$age_exit_chd-mULSAMFC$age_entry_chd
mULSAMFC$surv_chd20[mULSAMFC$surv_chd20>20] <- 20
mULSAMFC$surv_chd20[mULSAMFC$surv_chd20<0] <- 0.0001


## Make ratios 

mULSAMFC$LPC16_0.PC34_1 <-  log2(2^as.numeric(mULSAMFC$LPC18_0)/2^as.numeric(mULSAMFC$PC34_1))
mULSAMFC$LPC16_0.PC34_2 <-  log2(2^as.numeric(mULSAMFC$LPC18_0)/2^as.numeric(mULSAMFC$PC34_2))

mULSAMFC$LPC18_1.PC34_1 <-  log2(2^as.numeric(mULSAMFC$LPC18_1)/2^as.numeric(mULSAMFC$PC34_1))
mULSAMFC$LPC18_1.PC36_1 <-  log2(2^as.numeric(mULSAMFC$LPC18_1)/2^as.numeric(mULSAMFC$PC36_1))

mULSAMFC$LPC18_2.PC34_2 <-  log2(2^as.numeric(mULSAMFC$LPC18_2)/2^as.numeric(mULSAMFC$PC34_2))
mULSAMFC$LPC18_2.PC36_2 <-  log2(2^as.numeric(mULSAMFC$LPC18_2)/2^as.numeric(mULSAMFC$PC36_2))

mULSAMFC$LPC18_0.PC36_1 <-  log2(2^as.numeric(mULSAMFC$LPC18_0)/2^as.numeric(mULSAMFC$PC36_1))
mULSAMFC$LPC18_0.PC36_2 <-  log2(2^as.numeric(mULSAMFC$LPC18_0)/2^as.numeric(mULSAMFC$PC36_2))

mULSAMFC$LPC18_2.hdl <-  log2(2^as.numeric(mULSAMFC$LPC18_2)/(mULSAMFC$hdl))
mULSAMFC$LPC18_2.ldl <-  log2(2^as.numeric(mULSAMFC$LPC18_2)/(mULSAMFC$ldl))
mULSAMFC$LPC18_1.hdl <-  log2(2^as.numeric(mULSAMFC$LPC18_1)/(mULSAMFC$hdl))
mULSAMFC$LPC18_1.ldl <-  log2(2^as.numeric(mULSAMFC$LPC18_1)/(mULSAMFC$ldl))


mark <- c("LPC16_0","LPC18_0","LPC18_1","LPC20_4","LPC16_0.PC34_1","LPC16_0.PC34_2","LPC18_0.PC36_1","LPC18_0.PC36_2", "LPC18_1.PC34_1","LPC18_1.PC36_1","LPC18_2.PC34_2","LPC18_2.PC36_2","LPC18_1.hdl","LPC18_1.ldl","LPC18_2.hdl","LPC18_2.ldl")



RESU <- NULL
RESU2 <- NULL
for (i in mark)
{
	mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(mULSAMFC[,i]))+sbp+bmi+smoke+antihyp+ldl+hdl+diab+log(tg), data =mULSAMFC)
	res <- c(i,paste(round(exp(as.numeric(mod$coefficients[2])),2)," (",round(exp(mod$coefficients[2]+qnorm(0.025)*summary(mod)$coefficients[2,3]),2),"-",round(exp(mod$coefficients[2]-qnorm(0.025)*summary(mod)$coefficients[2,3]),2),")",sep=""), 2*pnorm(-abs(as.numeric(mod$coefficients[2])/as.numeric(summary(mod)$coefficients[2,3]))))
	res2 <- c(i,as.numeric(mod$coefficients[2]), summary(mod)$coefficients[2,3])		
	RESU <- rbind(RESU,res)
	RESU2 <- rbind(RESU2,res2)
}






##################
#### TwinGene ####
##################

### Recode diabetes ###
mTWINGENEF$diab_baseline_new <- ifelse((mTWINGENEF$diabetess==1 & mTWINGENEF$diabetess_doctor==1) | mTWINGENEF$diabetess_medicine==1 | (mTWINGENEF$glucos>= 7 & mTWINGENEF$fasting==0) | (mTWINGENEF$glucos>= 11 & mTWINGENEF$fasting!=0) | (mTWINGENEF$diab_regd <= mTWINGENEF$check_d & !is.na(mTWINGENEF$diab_regd)),1,0)

mTWINGENEFchd <- mTWINGENEF[(mTWINGENEF$sub==1  | mTWINGENEF$incchd==1) & mTWINGENEF$fasting!=2,]



# CHD-free subjects
metabo_sub_nopchdT <- mTWINGENEFchd[mTWINGENEFchd$chd==0 | (mTWINGENEFchd$chd==1 & mTWINGENEFchd$chdd > mTWINGENEFchd$check_d),]



metabo_sub_nopchdT$LPC18_1.PC34_1 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_1)/2^as.numeric(metabo_sub_nopchdT$PC34_1))
metabo_sub_nopchdT$LPC18_1.PC36_1 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_1)/2^as.numeric(metabo_sub_nopchdT$PC36_1))


metabo_sub_nopchdT$LPC18_0.PC36_1 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_0)/2^as.numeric(metabo_sub_nopchdT$PC36_1))
metabo_sub_nopchdT$LPC18_0.PC36_2 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_0)/2^as.numeric(metabo_sub_nopchdT$PC36_2))

metabo_sub_nopchdT$LPC18_2.PC34_2 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_2)/2^as.numeric(metabo_sub_nopchdT$PC34_2))
metabo_sub_nopchdT$LPC18_2.PC36_2 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_2)/2^as.numeric(metabo_sub_nopchdT$PC36_2))


metabo_sub_nopchdT$LPC16_0.PC34_1 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_0)/2^as.numeric(metabo_sub_nopchdT$PC34_1))
metabo_sub_nopchdT$LPC16_0.PC34_2 <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_0)/2^as.numeric(metabo_sub_nopchdT$PC34_2))


metabo_sub_nopchdT$LPC18_2.hdl <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_2)/(metabo_sub_nopchdT$hdl))
metabo_sub_nopchdT$LPC18_2.ldl <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_2)/(metabo_sub_nopchdT$ldl))
metabo_sub_nopchdT$LPC18_1.hdl <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_1)/(metabo_sub_nopchdT$hdl))
metabo_sub_nopchdT$LPC18_1.ldl <-  log2(2^as.numeric(metabo_sub_nopchdT$LPC18_1)/(metabo_sub_nopchdT$ldl))



library(Hmisc)
markT <- c("LPC16_0","LPC18_0","LPC18_1","LPC20_4","LPC16_0.PC34_1","LPC16_0.PC34_2","LPC18_0.PC36_1","LPC18_0.PC36_2", "LPC18_1.PC34_1","LPC18_1.PC36_1","LPC18_2.PC34_2","LPC18_2.PC36_2","LPC18_1.hdl","LPC18_1.ldl","LPC18_2.hdl","LPC18_2.ldl")



metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$tg) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi) & !is.na(metabo_sub_nopchdT$ldl),]


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
MAIN2 <- NULL
STRATAM <- NULL
for (i in markT)
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
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
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
		tg=mean(log(metabo_sub_nopchdT4$tg)),
		ldl=mean(metabo_sub_nopchdT4$ldl),
		hdl=mean(metabo_sub_nopchdT4$hdl),
		diab_baseline_new=mean(metabo_sub_nopchdT4$diab_baseline_new))

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
#### SO use STRATAM only for M520.340T346.671 ####
REST2 <- rbind(MAIN2[1,],MAIN2[2,],STRATAM[3,],MAIN2[4,],MAIN2[5,],MAIN2[6,],MAIN2[7,],MAIN2[8,],MAIN2[9,],MAIN2[10,],MAIN2[11,],MAIN2[12,],STRATAM[13,],MAIN2[14,],STRATAM[15,],MAIN2[16,])



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

chd_meta <- cbind(markT,paste(round(exp(meta_b),2)," (",round(exp(meta_b+qnorm(0.025)*meta_s),2),"-",round(exp(meta_b-qnorm(0.025)*meta_s),2),")",sep=""),2*pnorm(-abs(meta_b/meta_s)))


