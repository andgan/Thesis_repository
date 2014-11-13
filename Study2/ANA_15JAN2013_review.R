
#----------------------------------------------
# Filename: Analysis_15JAN2013.txt
# Study: Genscore
# Author: Andrea Ganna
# Date: 23FEB2013
# Updated: 17APR2013 - We added extra analysis for review, 2nd round.
#          24APR2013 - Added the cardiogrampluscd4 polygene score (~80000) SNPs, but this analysis is not included in the final paper
# Purpose: Analysis of the data for the genscore paper, after review
# Note: This are the changes done to answer review issues:
# 1. Check proportionality hazard
# 2. Add a new MGRS score, without CHD-associated SNPs
# 3. Check the biased-corrected clinical NRI
# 4. Check the variance explained by POLYGENE and GWAS-significant
# 5. Check risk factor distribution in individuals without risk factors and with at least one missing risk factor.
# 6. Calculate the NRI with the new risk categories <10,10-20,>20 %
# 7. Modify Figure 2 with new risk categories
#
# SECOND REVIEW:
# - New changes: add a forest plot
#-----------------------------------------------
# Data used: ULSAM_GOSH_TWGE_final_review.csv single_ana.txt single_ana_chd.txt single_ana_info.txt annotation.txt NetBen.r
# Data created: suppl_fig_1.jpeg single_ana.csv figure_2.jpeg suppl_figure_2a.jpeg suppl_figure_2b.jpeg suppl_figure_3.jpeg
#-----------------------------------------------
# OP: R 2.12.1
#-----------------------------------------------*/

####### SUMMARY ########

### 1. Import function and prepare data
### 2. DESCRIPTIVE and ASSUMPTIONS
###    a. TEST PROPORTIONALITY HAZARD
###    b. BASELINE CHARACTERISTICS
### 3. POLYGENE SCORE (suppl_fig_1.jpeg)
### 4. ASSOCIATION
###    a. CHD 
###    b. SINGLE SNP ASSOCIATION (single_ana.csv)
###    c. FHS risk factors 
### 5. PREDICTION MEASURES 
###    a. ALL-MEASURES
###    b. RECLASSIFICATION TABLE TO CALCULATE RECLASSIFICATION IN EVENTS AND NON-EVENTS
###    c. CALCULATE PREVENTION AMONG SCREENED
### 6. FIGURES
###    a. Figure 2 (figure_2.jpeg)
###    b. SUPPLEMENTARY FIGURE 2 (suppl_figure_2a.jpeg, suppl_figure_2b.jpeg)
###    c. SUPPLEMENTARY FIGURE 3 (suppl_figure_3.jpeg)
### 7. OTHER ANALYSIS
###    a. C-STATISTIC FROM LOGISTIC
###    b. NET BENEFIT
### SECOND REVIEW
### Forest plot
### Check association of 100 random scores

############################################
############################################
### 1. Import function and prepare data ####
############################################
############################################


#LIBRARY
library(rms)
library(survival)
library(gdata)
library(meta)
library(Hmisc)
library(aod)
library(survC1)


## LOAD FUNCTIONS

#NRI CLASSIC
NRI_classic=function(predwith,predwithout,out,cat){
P1 <- predwith
P2 <- predwithout
cutP1_e <-cut(P1[out==1], breaks=cat)
cutP1_ne <-cut(P1[out==0], breaks=cat)
cutP2_e <-cut(P2[out==1], breaks=cat)
cutP2_ne <-cut(P2[out==0], breaks=cat)
comb_e <- table (cutP2_e, cutP1_e)
comb_ne <- table (cutP2_ne, cutP1_ne)
e<- sum(comb_e)
ne <- sum(comb_ne)
up_e<- sum(upperTriangle(comb_e))/e
down_e<- sum(lowerTriangle(comb_e))/e
up_ne<- sum(upperTriangle(comb_ne))/ne
down_ne<- sum(lowerTriangle(comb_ne))/ne
NRI <- ((up_e-down_e)-(up_ne-down_ne))*100
seNRI <- ((up_e+down_e)/e+(up_ne+down_ne)/ne)^0.5
zNRI <- (NRI/100)/seNRI
pNRI <- round(2 * (pnorm(-abs(zNRI))),4)
return(data.frame(NRI,seNRI,pNRI))}

#clinicalNRI
NRI_clinical=function(predwith,predwithout,out,cat){
P1 <- predwith
P2 <- predwithout
cutP1_e <-cut(P1[out==1], breaks=cat)
cutP1_ne <-cut(P1[out==0], breaks=cat)
cutP2_e <-cut(P2[out==1], breaks=cat)
cutP2_ne <-cut(P2[out==0], breaks=cat)
comb_e <- table (cutP2_e, cutP1_e)
comb_ne <- table (cutP2_ne, cutP1_ne)
e<- sum(comb_e[2,])
ne <- sum(comb_ne[2,])
up_e<- comb_e[2,3]/e
down_e<- comb_e[2,1]/e
up_ne<- comb_ne[2,3]/ne
down_ne<- comb_ne[2,1]/ne
NRI <- ((up_e-down_e)-(up_ne-down_ne))*100
seNRI <- ((up_e+down_e)/e+(up_ne+down_ne)/ne)^0.5
zNRI <- (NRI/100)/seNRI
pNRI <- round(2 * (pnorm(-abs(zNRI))),4)
return(data.frame(NRI,seNRI,pNRI))}

#IDI
IDI=function(predwith,predwithout,out){
P1 <- predwith
P2 <- predwithout
IDI <- (mean(P1[out==1])-mean(P2[out==1]))-(mean(P1[out==0])-mean(P2[out==0]))
d <- P1-P2
seIDI <- sqrt((sd(P1[out==1] - P2[out==1])^2 / length(d[out==1])) + (sd(P1[out==0] - P2[out==0])^2 / length (d[out==0])))
zIDI <- IDI/seIDI
pIDI <- round(2 * (pnorm(-abs(zIDI))),4)
return(data.frame(IDI,seIDI,pIDI))}


# IMPORT data
setwd("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results")

d <- read.csv("ULSAM_GOSH_TWGE_final_review.csv")

# CHECK MISSING PATTERN
d$check_d <- as.Date(d$check_d,"%d%B%Y")
d_ <- d[!is.na(d$check_d),]

fin <- NULL
for (i in c("satsa","gender","harmony","octo","twingene")){
a <- sum(is.na(d_$bmi[d_$study==i]))/length(d_$bmi[d_$study==i])
b <-sum(is.na(d_$sbp[d_$study==i]))/length(d_$sbp[d_$study==i])
c <-sum(is.na(d_$hdl[d_$study==i]))/length(d_$hdl[d_$study==i])
e <-sum(is.na(d_$tc[d_$study==i]))/length(d_$tc[d_$study==i])
f <-sum(is.na(d_$antihyp[d_$study==i]))/length(d_$antihyp[d_$study==i])
g <-sum(is.na(d_$smoke[d_$study==i]))/length(d_$smoke[d_$study==i])
h <-sum(is.na(d_$diab[d_$study==i]))/length(d_$diab[d_$study==i])
fin <- rbind(fin,c(i,a,b,c,e,f,g,h))
}

# Dates

d$check_d <- as.Date(d$check_d,"%d%B%Y")
d$chdd <- as.Date(d$chdd,"%d%B%Y")
d$cvdd <- as.Date(d$cvdd,"%d%B%Y")
d$isd <- as.Date(d$isd,"%d%B%Y")
d$hsd <- as.Date(d$hsd,"%d%B%Y")
d$hfd <- as.Date(d$hfd,"%d%B%Y")
d$cadd <- as.Date(d$cadd,"%d%B%Y")

## Determine incident data and study_specific

d_inc_chd <- d[d$chdd > d$check_d & !is.na(d$age_entry_chd) & !is.na(d$SEX) & !is.na(d$tc) & !is.na(d$hdl) & !is.na(d$smoke) & !is.na(d$antihyp) & !is.na(d$sbp) & !is.na(d$diab) & !is.na(d$bmi),]

d_inc_cad <- d[d$cadd > d$check_d & !is.na(d$age_entry_chd) & !is.na(d$SEX) & !is.na(d$tc) & !is.na(d$hdl) & !is.na(d$smoke) & !is.na(d$antihyp) & !is.na(d$sbp) & !is.na(d$diab) & !is.na(d$bmi),]

d_inc_is <- d[d$isd > d$check_d & !is.na(d$SEX),]
d_inc_hs <- d[d$hsd > d$check_d & !is.na(d$SEX),]
d_inc_hf <- d[d$hfd > d$check_d & !is.na(d$SEX),]

### ADD POLYGENE SCORE (CARDIOGRAM, NOT USED IN FINAL PAPER)####
tw <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_twingene_cardio4cds.txt", stringsAsFactor=F)
ul <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_ulsam_cardio4cds.txt", stringsAsFactor=F)
ul$V18 <- paste("UL",ul$V18, sep="")
go <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_gosh_cardio4cds.txt", stringsAsFactor=F)
fin <- rbind(tw,go,ul)
d_inc_chd <- merge(d_inc_chd,fin,by.x="GWAS_ID", by.y="V18", all.x=T)



#### CHECK IF DISTRIBUTION RISK FACTOR WAS DIFFERENT BEFORE INCLUDED AND EXCLUDED INDIVIDUALS ###

d_inc_chd_ne <- d[d$chdd > d$check_d & !is.na(d$age_entry_chd) & (is.na(d$SEX) | is.na(d$tc) & is.na(d$hdl) | is.na(d$smoke) | is.na(d$antihyp) | is.na(d$sbp) | is.na(d$diab) | is.na(d$bmi)),]

riskfmen<-aggregate(data.frame(d_inc_chd_ne$age,d_inc_chd_ne$tc,d_inc_chd_ne$sbp,d_inc_chd_ne$hdl,d_inc_chd_ne$bmi),by=list(d_inc_chd_ne$SEX),mean, na.rm=T)
riskfsd<-aggregate(data.frame(d_inc_chd_ne$age,d_inc_chd_ne$tc,d_inc_chd_ne$sbp,d_inc_chd_ne$hdl,d_inc_chd_ne$bmi),by=list(d_inc_chd_ne$SEX),sd, na.rm=T)
riskcon <- xtabs(cbind(smoke,antihyp,diab,incchd)~SEX, data=d_inc_chd)

nrow(d_inc_chd_ne[d_inc_chd_ne$SEX==1,])
nrow(d_inc_chd_ne[d_inc_chd_ne$SEX==2,])


table(d_inc_chd_ne$diab[d_inc_chd_ne$SEX==1])
table(d_inc_chd_ne$diab[d_inc_chd_ne$SEX==2])

table(d_inc_chd_ne$smoke[d_inc_chd_ne$SEX==1])
table(d_inc_chd_ne$smoke[d_inc_chd_ne$SEX==2])

table(d_inc_chd_ne$antihyp[d_inc_chd_ne$SEX==1])
table(d_inc_chd_ne$antihyp[d_inc_chd_ne$SEX==2])


t.test(d_inc_chd_ne$age[d_inc_chd_ne$SEX==1],d_inc_chd$age[d_inc_chd$SEX==1])
t.test(d_inc_chd_ne$age[d_inc_chd_ne$SEX==2],d_inc_chd$age[d_inc_chd$SEX==2])

t.test(d_inc_chd_ne$tc[d_inc_chd_ne$SEX==1],d_inc_chd$tc[d_inc_chd$SEX==1])
t.test(d_inc_chd_ne$tc[d_inc_chd_ne$SEX==2],d_inc_chd$tc[d_inc_chd$SEX==2])

t.test(d_inc_chd_ne$sbp[d_inc_chd_ne$SEX==1],d_inc_chd$sbp[d_inc_chd$SEX==1])
t.test(d_inc_chd_ne$sbp[d_inc_chd_ne$SEX==2],d_inc_chd$sbp[d_inc_chd$SEX==2])

t.test(d_inc_chd_ne$hdl[d_inc_chd_ne$SEX==1],d_inc_chd$hdl[d_inc_chd$SEX==1])
t.test(d_inc_chd_ne$hdl[d_inc_chd_ne$SEX==2],d_inc_chd$hdl[d_inc_chd$SEX==2])

t.test(d_inc_chd_ne$bmi[d_inc_chd_ne$SEX==1],d_inc_chd$bmi[d_inc_chd$SEX==1])
t.test(d_inc_chd_ne$bmi[d_inc_chd_ne$SEX==2],d_inc_chd$bmi[d_inc_chd$SEX==2])

chisq.test(t(rbind(table(d_inc_chd_ne$smoke[d_inc_chd_ne$SEX==1]),table(d_inc_chd$smoke[d_inc_chd$SEX==1]))))
chisq.test(t(rbind(table(d_inc_chd_ne$smoke[d_inc_chd_ne$SEX==2]),table(d_inc_chd$smoke[d_inc_chd$SEX==2]))))

chisq.test(t(rbind(table(d_inc_chd_ne$diab[d_inc_chd_ne$SEX==1]),table(d_inc_chd$diab[d_inc_chd$SEX==1]))))
chisq.test(t(rbind(table(d_inc_chd_ne$diab[d_inc_chd_ne$SEX==2]),table(d_inc_chd$diab[d_inc_chd$SEX==2]))))


chisq.test(t(rbind(table(d_inc_chd_ne$antihyp[d_inc_chd_ne$SEX==1]),table(d_inc_chd$antihyp[d_inc_chd$SEX==1]))))
chisq.test(t(rbind(table(d_inc_chd_ne$antihyp[d_inc_chd_ne$SEX==2]),table(d_inc_chd$antihyp[d_inc_chd$SEX==2]))))

chisq.test(t(rbind(table(d_inc_chd_ne$antihyp),table(d_inc_chd$antihyp))))


d_inc_chd_ne <- d[d$chdd > d$check_d & !is.na(d$age_entry_chd) & (is.na(d$SEX) | is.na(d$tc) & is.na(d$hdl) | is.na(d$smoke) | is.na(d$antihyp) | is.na(d$sbp) | is.na(d$diab) | is.na(d$bmi)),]

# Check association between STAR_CHD and CHD in this subgroup with individuals with at least one missing risk factor
mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp , data=d_inc_chd_ne)

mod1 <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)



## Calculate individual risk at 10 years

d_inc_chd$survi <-  d_inc_chd$age_exit_chd - d_inc_chd$age_entry_chd


mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + bmi + study, data=d_inc_chd, surv=T, x=T, y=T)

mod2_lp <- predict.cph(mod2,type="lp")
d_inc_chd$Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

## For Var11 we also calculate quintiles

d_inc_chd$VAR11_q2 <- cut(d_inc_chd$VAR11,quantile(d_inc_chd$VAR11,(0:5)/5), include.lowest=TRUE, label=c(1,2,3,4,5))


## QUINTILES OF RISK
d_inc_chd$ALLCAT_weights_cl_wtccc_q2 <- cut(d_inc_chd$ALLCAT_weights_cl_wtccc,quantile(d_inc_chd$ALLCAT_weights_cl_wtccc,(0:5)/5), include.lowest=TRUE, label=c(1,2,3,4,5))
d_inc_chd$ALLCAT_chd_weights_cl_wtccc_q2 <- cut(d_inc_chd$ALLCAT_chd_weights_cl_wtccc,quantile(d_inc_chd$ALLCAT_chd_weights_cl_wtccc,(0:5)/5), include.lowest=TRUE, label=c(1,2,3,4,5))
d_inc_chd$VAR9_q <- cut(d_inc_chd$VAR9,quantile(d_inc_chd$VAR9,(0:4)/4), include.lowest=TRUE, label=c(1,2,3,4))
d_inc_chd$V17_q <- cut(d_inc_chd$V17,quantile(d_inc_chd$VAR9,(0:4)/4), include.lowest=TRUE, label=c(1,2,3,4))


# Individual studies

d_inc_chd_gosh <- d_inc_chd[d_inc_chd$study %in% c("gender","harmony","octo","satsa"),]
d_inc_chd_ulsam <- d_inc_chd[d_inc_chd$study=="ulsam",]
d_inc_chd_twingene <- d_inc_chd[d_inc_chd$study=="twingene",]


## Add another score FHS+CHD_star
d_inc_chd$fhs_star <- d_inc_chd$ALLCAT_fhs_cl_wtccc+d_inc_chd$STAR_chd


## Variable lists

l_sens <- c("ALLCAT_chd","ALLCAT","STAR_chd")
l_sens_q <- c("ALLCAT_chd_q","ALLCAT_q","STAR_chd_q")


l_s <- c("ALLCAT_weights_cl_wtccc","ALLCAT_chd_weights_cl_wtccc","VAR9","ALLCAT_bmi_cl_wtccc","ALLCAT_hdl_cl_wtccc","ALLCAT_sbp_cl_wtccc","ALLCAT_smoke_cl_wtccc","ALLCAT_t2d_cl_wtccc","ALLCAT_tc_cl_wtccc","ALLCAT_fhs_cl_wtccc","STAR_chd","fhs_star")

l_q <- c("ALLCAT_weights_cl_wtccc_q","ALLCAT_chd_weights_cl_wtccc_q","V17_q","STAR_chd_q")

l_p <- c("ALLCAT_weights_cl_wtccc","ALLCAT_chd_weights_cl_wtccc","V17")


#######################################
#######################################
### 2. DESCRIPTIVE and ASSUMPTIONS ####
#######################################
#######################################

#######################################
###  a. TEST PROPORTIONALITY HAZARD ###
#######################################

mod_overall<- coxph (Surv(age_entry_chd,age_exit_chd, incchd==1) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + bmi +study, data=d_inc_chd)
zph_overall <- cox.zph(mod_overall)
plot(zph_overall[2])

## With survival
mod_overall<- coxph (Surv(survi, incchd==1) ~  age + SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + bmi +study, data=d_inc_chd)
zph_overall <- cox.zph(mod_overall)
plot(zph_overall[3])



mod_gosh <- coxph (Surv(age_entry_chd,age_exit_chd, incchd==1) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + bmi, data=d_inc_chd_gosh)
zph_gosh <- cox.zph(mod_gosh)
plot(zph_gosh[2])

mod_ulsam <- coxph (Surv(age_entry_chd,age_exit_chd, incchd==1) ~  STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + bmi, data=d_inc_chd_ulsam)
zph_ulsam <- cox.zph(mod_ulsam)
plot(zph_ulsam[2])


mod_twingene <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd==1) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + bmi, data=d_inc_chd_twingene),d_inc_chd_twingene$GWAS_FID)
zph_twingene <- cox.zph(mod_twingene)
plot(zph_twingene[2])


zph_twingene <- cox.zph(mod_twingene)
plot(zph_twingene[2])
zph_twingene <- cox.zph(mod_twingene, transform="rank")
plot(zph_twingene[2])
zph_twingene <- cox.zph(mod_twingene, transform="log")
plot(zph_twingene[2])
zph_twingene <- cox.zph(mod_twingene, transform="identity")
plot(zph_twingene[2])


####################################
###  b. BASELINE CHARACTERISTICS ###
####################################


nrow(d_inc_chd)
median(d_inc_chd$survi)
quantile(d_inc_chd$survi,(0:4)/4)
table(d_inc_chd$study)

table(d_inc_chd$incchd[d_inc_chd$SEX==1])
table(d_inc_chd$incchd[d_inc_chd$SEX==2])

# Risk factors distribution
riskfmen<-aggregate(data.frame(d_inc_chd$age,d_inc_chd$tc,d_inc_chd$sbp,d_inc_chd$hdl,d_inc_chd$bmi),by=list(d_inc_chd$SEX),mean, na.rm=T)
riskfsd<-aggregate(data.frame(d_inc_chd$age,d_inc_chd$tc,d_inc_chd$sbp,d_inc_chd$hdl,d_inc_chd$bmi),by=list(d_inc_chd$SEX),sd, na.rm=T)
riskcon <- xtabs(cbind(smoke,antihyp,diab,incchd)~SEX, data=d_inc_chd)

nrow(d_inc_chd[d_inc_chd$SEX==1,])
nrow(d_inc_chd[d_inc_chd$SEX==2,])


#Association
ass <-robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab + bmi + study, data=d_inc_chd),d_inc_chd$GWAS_FID)


median(d_inc_chd_twingene$survi)
quantile(d_inc_chd_twingene$survi,(0:4)/4)

tt <- rbind(d_inc_chd_ulsam, d_inc_chd_gosh)

median(tt$survi)
quantile(tt$survi,(0:4)/4)
table(tt$incchd)



##########################
##########################
### 3. POLYGENE SCORE ###
#########################
#########################

# To understand which is the best cut_off to use in the next analysis

# OVERALL


prof <- c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10","VAR11","VAR12","VAR13","VAR14","VAR15","VAR16","VAR17")

R2_e <- NULL
C_e <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + SEX + scale(d_inc_chd[,i]) + study, data=d_inc_chd)
         bmod <- lrm(incchd ~  age + SEX  + study, data=d_inc_chd)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e <- c(R2_e,temp)
         C_e <-c(C_e,temp2)
         p_val<-c(p_val,anova(mod)["d_inc_chd","P"])}



# TWINGENE


prof <- c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10","VAR11","VAR12","VAR13","VAR14","VAR15","VAR16","VAR17")

R2_e_tw <- NULL
C_e_tw <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + SEX + scale(d_inc_chd_twingene[,i]) , data=d_inc_chd_twingene)
         bmod <- lrm(incchd ~  age + SEX, data=d_inc_chd_twingene)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_tw <- c(R2_e_tw,temp)
         C_e_tw <-c(C_e_tw,temp2)
         p_val<-c(p_val,anova(mod)["d_inc_chd_twingene","P"])}


# ULSAM 

d_inc_chd_u <- d_inc_chd[d_inc_chd$study %in% c("ulsam"),]

prof <- c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10","VAR11","VAR12","VAR13","VAR14","VAR15","VAR16","VAR17")

R2_e_u <- NULL
C_e_u <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + scale(d_inc_chd_u[,i]), data=d_inc_chd_u)
         bmod <- lrm(incchd ~   age , data=d_inc_chd_u)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_u <- c(R2_e_u,temp)
         C_e_u <-c(C_e_u,temp2)
         p_val<-c(p_val,anova(mod)["d_inc_chd_u","P"])}






# GOSH
d_inc_chd_g <- d_inc_chd[d_inc_chd$study %in% c("gender","harmony","octo","satsa"),]
prof <- c("VAR1","VAR2","VAR3","VAR4","VAR5","VAR6","VAR7","VAR8","VAR9","VAR10","VAR11","VAR12","VAR13","VAR14","VAR15","VAR16","VAR17")

R2_e_g <- NULL
C_e_g <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + scale(d_inc_chd_g[,i]) + SEX , data=d_inc_chd_g)
         bmod <- lrm(incchd ~   age + SEX , data=d_inc_chd_g)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_g <- c(R2_e_g,temp)
         C_e_g <-c(C_e_g,temp2)
         p_val<-c(p_val,anova(mod)["d_inc_chd_g","P"])}
         
         

jpeg("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/Paper/suppl_fig_1.jpeg",width=800,height=800,quality=100)


pvalues <- c(0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

plot(pvalues,C_e,xlab="P-values",ylab="C-statistic increment",type="l",lwd=2,lty=2,xaxt="n",ylim=c(-0.002,0.002))

#par(new=T)

#plot(pvalues,C_e_tw,xlab="",ylab="",type="l",lwd=1,xaxt="n",col="blue",ylim=c(0,0.016))

#par(new=T)

#plot(pvalues,C_e_g,xlab="",ylab="",type="l",lwd=1,xaxt="n",col="red",ylim=c(0,0.016))

#par(new=T)

#plot(pvalues,C_e_u,xlab="",ylab="",type="l",lwd=1,xaxt="n",col="green",ylim=c(0,0.016))



axis(side=1, at=c(0.00000005,0.0000005,0.000005,0.00005,0.0005,0.005,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),labels=c("5E-8",NA,NA,NA,NA,NA,"0.05","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"), xlim=c(0,1),las=2, tck=-.01)


#legend("topleft", inset=.05, col=c("black","blue","red","green"), lty=c(1,1,1,1), lwd=c(2,1,1,1), c("Overall", "llumina HumanOmniExpress ","Illumina Metabochip","Illumina Metabochip"), cex=0.8)

title("Polygene Score \n C-statistic Improvement Over Age, Sex and Study of Origin")

dev.off()

########### POLYGENE CARDIOGRAM  (NOT USED IN FINAL PAPER)  ########

tw <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_twingene_cardio4cds.txt", stringsAsFactor=F)

twm <- merge(d_inc_chd,tw,by.x="GWAS_ID", by.y="V18")

prof <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17")

R2_e_tw <- NULL
C_e_tw <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + SEX + scale(twm[,i]) , data=twm)
         bmod <- lrm(incchd ~  age + SEX, data=twm)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_tw <- c(R2_e_tw,temp)
         C_e_tw <-c(C_e_tw,temp2)
         p_val<-c(p_val,anova(mod)["twm","P"])}


ul <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_ulsam_cardio4cds.txt", stringsAsFactor=F)
ul$V18 <- paste("UL",ul$V18, sep="")

ulm <- merge(d_inc_chd,ul,by.x="GWAS_ID", by.y="V18")

prof <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17")

R2_e_ul <- NULL
C_e_ul <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age  + scale(ulm[,i]) , data=ulm)
         bmod <- lrm(incchd ~  age , data=ulm)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_ul <- c(R2_e_ul,temp)
         C_e_ul <-c(C_e_ul,temp2)
         p_val<-c(p_val,anova(mod)["ulm","P"])}
go <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/polygene/export_gosh_cardio4cds.txt", stringsAsFactor=F)

gom <- merge(d_inc_chd,go,by.x="GWAS_ID", by.y="V18")

prof <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17")

R2_e_go <- NULL
C_e_go <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + SEX + scale(gom[,i]) , data=gom)
         bmod <- lrm(incchd ~  age + SEX, data=gom)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_go <- c(R2_e_go,temp)
         C_e_go <-c(C_e_go,temp2)
         p_val<-c(p_val,anova(mod)["gom","P"])}


fin <- rbind(tw,go,ul)

finm <- merge(d_inc_chd,fin,by.x="GWAS_ID", by.y="V18")

prof <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17")

R2_e_fin <- NULL
C_e_fin <- NULL
p_val <- NULL
for (i in prof){
        mod <- lrm(incchd ~  age + SEX + scale(finm[,i]) +study, data=finm)
         bmod <- lrm(incchd ~  age + SEX + study, data=finm)
                        temp <- mod$stats[10] - bmod$stats[10]
                        temp2 <- mod$stats[6] - bmod$stats[6]
         R2_e_fin <- c(R2_e_fin,temp)
         C_e_fin <-c(C_e_fin,temp2)
         p_val<-c(p_val,anova(mod)["finm","P"])}




#######################
#######################
### 4. ASSOCIATION ####
#######################
#######################


##############
### a. CHD ###
##############


### Check association with this score without CHD, as asked by the reviewer ###
mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + allcat_no_chd_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)

### Check 20-years association in ULSAM
mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  STAR_chd + tc + hdl + smoke + antihyp + sbp + diab , data=d_inc_chd_ulsam)

cph(Surv(age_entry_chd,age_exit_chd, incchd) ~  STAR_chd + tc + hdl + smoke + antihyp + sbp + diab , data=d_inc_chd_ulsam)


### Check variance explained ####


## CHD-specific SNPs
mod1b <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)

mod1 <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)


((summary(mod1)$rsq[1]/summary(mod1)$rsq[2])-(summary(mod1b)$rsq[1]/summary(mod1b)$rsq[2]))*100

## Polygene score
mod1b <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)

mod1 <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + VAR9 + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)


((summary(mod1)$rsq[1]/summary(mod1)$rsq[2])-(summary(mod1b)$rsq[1]/summary(mod1b)$rsq[2]))*100



## Polygene score (GWAS Significant SNPs)
mod1b <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)

mod1 <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + VAR2 + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)


((summary(mod1)$rsq[1]/summary(mod1)$rsq[2])-(summary(mod1b)$rsq[1]/summary(mod1b)$rsq[2]))*100




### NORMAL ANALYSIS ###


# Continuous
res_c <- NULL
for (i in l_s){
	
mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + d_inc_chd[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_chd),d_inc_chd$GWAS_FID)
r <- c(i,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5])
res_c <- rbind(res_c,r)}


# Continuous (CAD)
res_c <- NULL
for (i in l_s){
	
mod <- robcov(coxph (Surv(age_entry_cad,age_exit_cad, inccad) ~  SEX + d_inc_cad[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=d_inc_cad),d_inc_cad$GWAS_FID)
r <- c(i,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5])
res_c <- rbind(res_c,r)}



# Ischeamic stroke
res_c <- NULL
for (i in l_s){
	
mod <- coxph (Surv(age_entry_is,age_exit_is, incis) ~  SEX + d_inc_is[,i] + study, data=d_inc_is)
r <- c(i,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5])
res_c <- rbind(res_c,r)}

# Hemmoragic stroke
res_c <- NULL
for (i in l_s){
	
mod <- coxph (Surv(age_entry_hs,age_exit_hs, inchs) ~  SEX + d_inc_hs[,i] + study, data=d_inc_hs)
r <- c(i,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5])
res_c <- rbind(res_c,r)}

# HF
res_c <- NULL
for (i in l_s){
	
mod <- coxph (Surv(age_entry_hf,age_exit_hf, inchf) ~  SEX + d_inc_hf[,i] + study, data=d_inc_hf)
r <- c(i,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5])
res_c <- rbind(res_c,r)}



# Quintiles and trend test

res_q <- NULL
for (i in l_q){
mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~ as.factor(d_inc_chd[,i])+  SEX  + tc + hdl + smoke + antihyp + sbp + diab + bmi + study, data=d_inc_chd),d_inc_chd$GWAS_FID)

zz <- c(1,2,3)
test.num <- zz %*% coef(mod)[1:3]
test.var <- zz %*% mod$var[1:3,1:3] %*% zz
z<-test.num/sqrt(test.var)
p <- 2*pnorm(-abs(z))

r <- c(i,summary(mod)$coefficients[3,2],exp(summary(mod)$coefficients[3,1]-qnorm(0.025)*summary(mod)$coefficients[3,3]),exp(summary(mod)$coefficients[3,1]+qnorm(0.025)*summary(mod)$coefficients[3,3]),summary(mod)$coefficients[3,5],p)
res_q <- rbind(res_q,r)}



#################################
### b. SINGLE SNP ASSOCIATION ###
#################################

single <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Data/From_plink/single_ana.txt", header=T)
single_chd <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Data/From_plink/single_ana_chd.txt", header=T)
single_info <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Data/From_plink/single_ana_info.txt", header=T)
anno <- read.table("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/SNPs_selection/annotation.txt", header=T, blank.lines.skip = FALSE)


s_def <- merge(d_inc_chd,single, by="GWAS_ID")
s_def_chd <- merge(d_inc_chd,single_chd, by="GWAS_ID")

res_tot <- NULL
for (i in single_info$SNPs)
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"BMI",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_tot <- rbind(res_tot,temp)}	

#write.table(res_tot,file="test.txt", quote=F, col.names=T, row.names=F, sep=",")


res_bmi <- NULL
for (i in single_info$SNPs[single_info$bmi==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"BMI",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_bmi <- rbind(res_bmi,temp)}	


# Chd are SNPs in CHD_star
res_chd <- NULL
for (i in ((length(s_def_chd)-45):length(s_def_chd)))
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def_chd[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def_chd),s_def_chd$GWAS_FID)
 freq <- ifelse(s_def_chd[,i] < 0.5,0,ifelse(s_def_chd[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(colnames(s_def_chd)[i],"_")[[1]][1],strsplit(colnames(s_def_chd)[i],"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"CHD",ifelse(colnames(s_def_chd)[i]%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(colnames(s_def_chd)[i]%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_chd <- rbind(res_chd,temp)}	

res_hdl <- NULL
for (i in single_info$SNPs[single_info$hdl==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"HDL",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_hdl <- rbind(res_hdl,temp)}	


res_sbp <- NULL
for (i in single_info$SNPs[single_info$sbp==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"SBP",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_sbp <- rbind(res_sbp,temp)}	


res_smoke <- NULL
for (i in single_info$SNPs[single_info$smoke==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"Smoke",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_smoke <- rbind(res_smoke,temp)}	


res_t2d <- NULL
for (i in single_info$SNPs[single_info$t2d==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"T2D",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_t2d <- rbind(res_t2d,temp)}	


res_tc <- NULL
for (i in single_info$SNPs[single_info$tc==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"TC",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_tc <- rbind(res_tc,temp)}	


res_fhs <- NULL
for (i in single_info$SNPs[single_info$fhs==1 & single_info$tc==0 & single_info$t2d==0 & single_info$smoke==0
& single_info$sbp==0 & single_info$hdl==0 & single_info$bmi==0])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"-",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_fhs <- rbind(res_fhs,temp)}	


res_allcat <- NULL
for (i in single_info$SNPs[single_info$fhs==0 & single_info$tc==0 & single_info$t2d==0 & single_info$smoke==0
& single_info$sbp==0 & single_info$hdl==0 & single_info$chd==0 & single_info$bmi==0 & single_info$allcat==1])
{mod <- robcov(coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + s_def[,i] + tc + hdl + smoke + antihyp + sbp + diab  + study, data=s_def),s_def$GWAS_FID)
 freq <- ifelse(s_def[,i] < 0.5,0,ifelse(s_def[,i] < 1.5,1,2))
 freq <- table(freq)
 alfreq <- (as.numeric(freq[3])+1/2*as.numeric(freq[2]))/sum(freq)		
 temp <- data.frame(strsplit(i,"_")[[1]][1],strsplit(i,"_")[[1]][2],alfreq,summary(mod)$coefficients[2,2],summary(mod)$coefficients[2,5],"-",ifelse(i%in%single_info$SNPs[single_info$fhs==1],1,0),ifelse(i%in%single_info$SNPs[single_info$allcat==1],1,0))
 colnames(temp) <- c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n")
 res_allcat <- rbind(res_allcat,temp)}	


res <- rbind(res_bmi,res_chd,res_hdl,res_sbp,res_smoke, res_t2d, res_tc, res_fhs, res_allcat)
fin <- merge(res, anno, by="SNPs")

final <- fin[,c("SNPs","Risk allele","Risk allele frequency","OR","P-value","Trait-specific MGRS","FHS MGRS y/n","Overall MGRS y/n","Disease_Trait","Reported_Gene_s_","PUBMEDID")]

write.csv(final,file="/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Results/single_ana.csv")

############################
### c. FHS risk factors  ###
############################

# Continuous 

### BMI ###
res_bmi <- NULL
for (i in l_s){
d_bmi <- d_inc_chd[!is.na(d_inc_chd$bmi),]	
mod <- robcov(ols (bmi ~  age + SEX + d_bmi[,i] + study,data=d_bmi, x=T, y=T), d_bmi$GWAS_FID)
r <- c(i,mod$coefficients[4],2*pnorm(-abs(mod$coefficients[4]/sqrt(diag(mod$var))[4])),mod$stats[4])
res_bmi <- rbind(res_bmi,r)}

mod_bmi_base <- robcov(ols (bmi ~  age + SEX + study,data=d_bmi, x=T, y=T), d_bmi$GWAS_FID)
R2_bmi <- (as.numeric(res_bmi[4,4]) - mod_bmi_base$stats[4])*100



### HDL ###

res_hdl <- NULL
for (i in l_s){
d_hdl <- d_inc_chd[!is.na(d_inc_chd$hdl),]	
mod <- robcov(ols (hdl ~  age + SEX + d_hdl[,i] + study + drug_lipid,data=d_hdl, x=T, y=T), d_hdl$GWAS_FID)
r <- c(i,mod$coefficients[4],2*pnorm(-abs(mod$coefficients[4]/sqrt(diag(mod$var))[4])),mod$stats[4])
res_hdl <- rbind(res_hdl,r)}

mod_hdl_base <- robcov(ols (hdl ~  age + SEX + study,data=d_hdl, x=T, y=T), d_hdl$GWAS_FID)
R2_hdl <- (as.numeric(res_hdl[5,4]) - mod_hdl_base$stats[4])*100


### SBP ###

res_sbp <- NULL
for (i in l_s){
d_sbp <- d_inc_chd[!is.na(d_inc_chd$sbp),]	
mod <- robcov(ols (sbp ~  age + SEX + d_sbp[,i] + study,data=d_sbp, x=T, y=T), d_sbp$GWAS_FID)
r <- c(i,mod$coefficients[4],anova(mod)["d_sbp","P"],mod$stats[4])
res_sbp <- rbind(res_sbp,r)}

mod_sbp_base <- robcov(ols (sbp ~  age + SEX + study,data=d_sbp, x=T, y=T), d_sbp$GWAS_FID)
R2_sbp <- (as.numeric(res_sbp[6,4]) - mod_sbp_base$stats[4])*100

### SMOKE ###

res_smoke <- NULL
for (i in l_s){
d_smoke <- d_inc_chd[!is.na(d_inc_chd$smoke),]	
mod <- robcov(lrm (smoke ~  age + SEX + d_smoke[,i] + study,data=d_smoke, x=T, y=T), d_smoke$GWAS_FID)
r <- c(i,mod$coefficients[4],anova(mod)["d_smoke","P"],mod$stats[10])
res_smoke <- rbind(res_smoke,r)}

mod_smoke_base <- robcov(lrm (smoke ~  age + SEX + study,data=d_smoke, x=T, y=T), d_smoke$GWAS_FID)
R2_smoke <- (as.numeric(res_smoke[7,4]) - mod_smoke_base$stats[10])*100


### DIABETES ###

res_diab <- NULL
for (i in l_s){
d_diab <- d_inc_chd[!is.na(d_inc_chd$diab),]	
mod <- robcov(lrm (diab ~  age + SEX + d_diab[,i] + study,data=d_diab, x=T, y=T), d_diab$GWAS_FID)
r <- c(i,mod$coefficients[4],anova(mod)["d_diab","P"],mod$stats[10])
res_diab <- rbind(res_diab,r)}

mod_diab_base <- robcov(lrm (diab ~  age + SEX + study,data=d_diab, x=T, y=T), d_diab$GWAS_FID)
R2_diab <- (as.numeric(res_diab[8,4]) - mod_diab_base$stats[10])*100


### TC ###
res_tc <- NULL
for (i in l_s){
d_tc <- d_inc_chd[!is.na(d_inc_chd$tc),]	
mod <- robcov(ols (tc ~  age + SEX + d_tc[,i] + study + drug_lipid,data=d_tc, x=T, y=T), d_tc$GWAS_FID)
r <- c(i,mod$coefficients[4],2*pnorm(-abs(mod$coefficients[4]/sqrt(diag(mod$var))[4])),mod$stats[4])
res_tc <- rbind(res_tc,r)}

mod_tc_base <- robcov(ols (tc ~  age + SEX + study,data=d_tc, x=T, y=T), d_tc$GWAS_FID)
R2_tc <- (as.numeric(res_tc[9,4]) - mod_tc_base$stats[4])*100



	
### INDIVIDUAL RISK ###
d_Rt_wo <- d_inc_chd[!is.na(d_inc_chd$Rt_wo) & d_inc_chd$study!="twingene",]
mod2 <- cph (Surv(survi, incchd) ~ SEX +  tc + hdl + smoke + antihyp + sbp + diab + age + study , data=d_Rt_wo, surv=T, x=T, y=T)
mod2_lp <- predict.cph(mod2,type="lp")
d_Rt_wo$Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

res_Rt_wo <- NULL
for (i in l_s){
mod <- robcov(ols (Rt_wo ~  age + SEX + d_Rt_wo[,i] + study,data=d_Rt_wo, x=T, y=T), d_Rt_wo$GWAS_FID)
r <- c(i,mod$coefficients[4],anova(mod)["d_Rt_wo","P"],mod$stats[4])
res_Rt_wo <- rbind(res_Rt_wo,r)}


###############################
###############################
### 5. PREDICTION MEASURES ####
###############################
###############################



#######################
### a. ALL-MEASURES ###
#######################

res_nri <- NULL
d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]


for (i in l_s){
				
## INCLUDING GENETIC SCORE

mod1 <- cph (Surv(survi, incchd) ~   SEX + scale(d_inc_chd_p[,i]) + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)

mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)

# GOF
cutRt <- as.factor(cut(log(Rt_w), quantile(log(Rt_w),(0:5)/5), labels=1:5, include.lowest=T))
mod1_gof <- cph (Surv(survi, incchd) ~  SEX + scale(d_inc_chd_p[,i]) + tc + hdl + smoke + antihyp + sbp + diab + age + study + cutRt, data=d_inc_chd_p, surv=T, x=T, y=T)

GOF<- wald.test(u <- vcov(mod1_gof),coef(mod1_gof), Terms=14:17)
GOF<- as.vector(GOF$result$chi2)[3]


### WITHOUT GENETIC SCORE

# Model
mod2 <- cph (Surv(survi, incchd) ~ SEX +  tc + hdl + smoke + antihyp + sbp + diab + age + study , data=d_inc_chd_p, surv=T, x=T, y=T)

mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)



# Study with 10 years of follow-up
d_inc_chd_10 <- d_inc_chd_p
d_inc_chd_10$incchd <- ifelse(d_inc_chd_p$survi>10,0,d_inc_chd_p$incchd)
d_inc_chd_10$survi <- ifelse(d_inc_chd_p$survi>10,10,d_inc_chd_p$survi)



### NRI ###
nri <- NRI_classic(Rt_w,Rt_wo,d_inc_chd_10$incchd,c(0,0.10,0.20,100))

### CLINICAL NRI ###
clnri <- NRI_clinical(Rt_w,Rt_wo,d_inc_chd_10$incchd,c(0,0.05,0.20,100))

### C-INDEX ###
c_w <- 1-rcorr.cens(Rt_w,Surv(d_inc_chd_10$survi, d_inc_chd_10$incchd))[1]
c_wo <- 1-rcorr.cens(Rt_wo,Surv(d_inc_chd_10$survi, d_inc_chd_10$incchd))[1]
c <- rcorrp.cens(Rt_w,Rt_wo,Surv(d_inc_chd_10$survi, d_inc_chd_10$incchd)) 

c_m <- c_w - c_wo
					
r <- c(i,nri[1],((nri[1]/100)+qnorm(0.025)*nri[2])*100,((nri[1]/100)-qnorm(0.025)*nri[2])*100, nri[3],clnri[1],((clnri[1]/100)+qnorm(0.025)*clnri[2])*100,((clnri[1]/100)-qnorm(0.025)*clnri[2])*100,clnri[3],c_m,GOF)
res_nri <- rbind(res_nri,r)
}

############################################################################################
#### b. RECLASSIFICATION TABLE TO CALCULATE RECLASSIFICATION IN EVENTS AND NON-EVENTS ####
############################################################################################


d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]

# Study with 10 years of follow-up
d_inc_chd_10 <- d_inc_chd_p
d_inc_chd_10$incchd <- ifelse(d_inc_chd_p$survi>10,0,d_inc_chd_p$incchd)
d_inc_chd_10$survi <- ifelse(d_inc_chd_p$survi>10,10,d_inc_chd_p$survi)



## TO BE DONE DIRECTELY IN THE TERMINAL< OTHERWISE DOES NOT WORK

## STAR_CHD
library(PredictABEL)
mod1 <- cph (Surv(survi, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)

mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

reclassification(data= d_inc_chd_10, cOutcome=27, predrisk1=Rt_wo, predrisk2= Rt_w, cutoff=c(0,.10,.20,1))

((84+63)/2626)-((109+92)/2626)
((10+20)/388)-((11+8)/388)


## OVERALL MGRS
mod1 <- cph (Surv(survi, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)

mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

reclassification(data= d_inc_chd_10, cOutcome=27, predrisk1=Rt_wo, predrisk2= Rt_w, cutoff=c(0,.10,.20,1))

((56+48)/2626)-((66+73)/2626)
((7+13)/388)-((5+4)/388)


## POLYGENE SCORE
mod1 <- cph (Surv(survi, incchd) ~  SEX + VAR9 + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)

mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)
mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

reclassification(data= d_inc_chd_10, cOutcome=27, predrisk1=Rt_wo, predrisk2= Rt_w, cutoff=c(0,.10,.20,1))

((16+15)/2626)-((20+11)/2626)
((2+5)/388)-((1+1)/388)


#### TRY THE CLINICAL NRI (SEE EXCEL SPREADSHEET IN THE Paper/ATVB/Review/Clinical_NRI)



###############################################
#### c. CALCULATE PREVENTION AMONG SCREENED ###
###############################################


d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]

d_inc_chd_p <- d_inc_chd


SS <- NULL
for (i in 1:100)
{
	set.seed(123+i)
	d_inc_chd_b <- d_inc_chd_p[sample(nrow(d_inc_chd_p), replace=T),]
	
	# Study with 10 years of follow-up
	d_inc_chd_10 <- d_inc_chd_b
	d_inc_chd_10$incchd <- ifelse(d_inc_chd_b$survi>10,0, d_inc_chd_b$incchd)
	d_inc_chd_10$survi <- ifelse(d_inc_chd_b$survi>10,10, d_inc_chd_b$survi)

	
	mod1 <- cph (Surv(survi, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + age + study, data= d_inc_chd_b, surv=T, x=T, y=T)
	mod1_lp <- predict.cph(mod1,type="lp")
	Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)

	mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data= d_inc_chd_b, surv=T, x=T, y=T)
	mod2_lp <- predict.cph(mod2,type="lp")
	Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)

	recl<-length(Rt_wo[Rt_wo>0.1 & Rt_wo<0.2])
	prev <- sum((Rt_wo>0.1 & Rt_wo<0.2) & (Rt_w>0.2) & d_inc_chd_10$incchd==1)*0.20
	ss <- recl/prev
	SS <- c(SS,ss)
	print(prev)
}

SS[is.infinite(SS)] <- NA

mean(SS,na.rm=T)+1.96*sd(SS,na.rm=T)
mean(SS,na.rm=T)-1.96*sd(SS,na.rm=T)


1272/4




###################
###################
### 6. FIGURES ####
###################
###################


#####################
#### a. FIGURE 2 ####
#####################

d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]


#### CHD-specific MGRS ####

mod1 <- cph (Surv(survi, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)

mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)


### WITHOUT GENETIC SCORE

# Model
mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)

mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)


# Study with 10 years of follow-up
d_inc_chd_10 <- d_inc_chd_p
d_inc_chd_10$incchd <- ifelse(d_inc_chd_p$survi>10,0,d_inc_chd_p$incchd)
d_inc_chd_10$survi <- ifelse(d_inc_chd_p$survi>10,10,d_inc_chd_p$survi)



cutRt_w <-as.numeric(cut(Rt_w, breaks=c(0,0.10,0.20,100),labels=c(1,2,3)))
cutRt_wo <-as.numeric(cut(Rt_wo, breaks=c(0,0.10,0.20,100),labels=c(1,2,3)))



col_ca <-ifelse(cutRt_w[d_inc_chd_10$incchd==1] > cutRt_wo[d_inc_chd_10$incchd==1],"blue",
ifelse(cutRt_w[d_inc_chd_10$incchd==1] < cutRt_wo[d_inc_chd_10$incchd==1],"red","black"))

col_co<-ifelse(cutRt_w[d_inc_chd_10$incchd==0] > cutRt_wo[d_inc_chd_10$incchd==0],"red",
ifelse(cutRt_w[d_inc_chd_10$incchd==0] < cutRt_wo[d_inc_chd_10$incchd==0],"blue","darkgrey"))



cutP1_e <-cut(Rt_w[d_inc_chd_10$incchd==1], breaks=c(0,0.10,0.20,100))
cutP1_ne <-cut(Rt_w[d_inc_chd_10$incchd==0], breaks=c(0,0.10,0.20,100))
cutP2_e <-cut(Rt_wo[d_inc_chd_10$incchd==1], breaks=c(0,0.10,0.20,100))
cutP2_ne <-cut(Rt_wo[d_inc_chd_10$incchd==0], breaks=c(0,0.10,0.20,100))
comb_e <- table (cutP2_e, cutP1_e)
comb_ne <- table (cutP2_ne, cutP1_ne)


jpeg("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/Paper/figure_2.jpeg",width=600,height=600,quality=100)


plot(Rt_wo[d_inc_chd_10$incchd==0]*100,Rt_w[d_inc_chd_10$incchd==0]*100, log="xy", pch=21, col=col_co,xaxt="n", yaxt="n", axes=F,xlab="",ylab="", xlim=c(1,100),ylim=c(1,100))
par(new=T)
plot(Rt_wo[d_inc_chd_10$incchd==1]*100,Rt_w[d_inc_chd_10$incchd==1]*100, log="xy", pch=4, col=col_ca,xaxt="n", yaxt="n", axes=F,
,xlab="",ylab="", xlim=c(1,100),ylim=c(1,100))
abline(v = 10)
abline(v = 20)
abline(h = 10)
abline(h = 20)
mtext(side = 1, text = "10-Year Risk Predicted by model including only FHS risk factors", line = 2.2 ,cex=0.85)
mtext(side = 2, text = "10-Year Risk Predicted by model including \n FHS risk factors and CHD-specific MGRS", line = 2.2, cex=0.85)
axis(1,at=c(1,3,5,10,20,40,70,100), las=1, label=c(1,3,5,10,20,40,70,100), col='black', lwd=1)
axis(2,at=c(1,3,5,10,20,40,70,100), las=1, label=c(01,3,5,10,20,40,70,100), col='black', lwd=1)

legend(1.5,70,legend=c("CHD event","Non-event"), col=c("black", "darkgrey"), pch=c(4,21), cex=c(0.8), bg="white")

#text(3,7,paste(comb_e[1,2],"/",comb_ne[1,2],sep=""), cex=0.8)
#text(7,2.5,paste(comb_e[2,1],"/",comb_ne[2,1],sep=""), cex=0.8)
#text(6.3,15,paste(comb_e[2,3],"/",comb_ne[2,3],sep=""), cex=0.8)
#text(15,6.5,paste(comb_e[3,2],"/",comb_ne[3,2],sep=""), cex=0.8)
#text(13,30,paste(comb_e[3,4],"/",comb_ne[3,4],sep=""), cex=0.8)
#text(30,13,paste(comb_e[4,3],"/",comb_ne[4,3],sep=""), cex=0.8)




dev.off()

###################################
#### b .SUPPLEMENTARY FIGURE 2 ####
###################################

d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]

mod1 <- cph (Surv(survi, incchd) ~  SEX + STAR_chd + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)


mod1_lp <- predict.cph(mod1,type="lp")
Rt_w <- (1-survest(mod1, linear.predictors=mod1_lp,times=10, conf.int=F)$surv)


### WITHOUT GENETIC SCORE

# Model
mod2 <- cph (Surv(survi, incchd) ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd_p, surv=T, x=T, y=T)

mod2_lp <- predict.cph(mod2,type="lp")
Rt_wo <- (1-survest(mod2, linear.predictors=mod2_lp,times=10, conf.int=F)$surv)


# Study with 10 years of follow-up
d_inc_chd_10 <- d_inc_chd_p
d_inc_chd_10$incchd <- ifelse(d_inc_chd_p$survi>10,0,d_inc_chd_p$incchd)
d_inc_chd_10$survi <- ifelse(d_inc_chd_p$survi>10,10,d_inc_chd_p$survi)


jpeg("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/Paper/suppl_figure_2a.jpeg",width=800,height=800,quality=100)


par(mar=c(4,5,3,3))
x <- Rt_wo[d_inc_chd_10$incchd==1]
y <- Rt_w[d_inc_chd_10$incchd==1]-Rt_wo[d_inc_chd_10$incchd==1]

r <- lowess(x[x<0.5],y[x<0.5],f=0.8)
plot(x,y,col="antiquewhite4",xlim=c(0,0.8),ylim=c(-0.10,+0.10),ylab="",xlab="",axes=F)
title("CHD events",xlab="10-year risk predicted by model containing only FHS risk factors",ylab="Differences in 10-year risk between model \n with FHS risk factors + overall MGRS and only FHS", cex.lab=0.85)
lines(r,lwd=3)
axis(1,at=c(0,0.05,0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80),label=c(0,5,10,15,20,30,40,50,60,70,80))
axis(2,at=c(-0.10,-0.05,0,0.05,0.10),label=c(-10,-5,0,+5,+10),las=1)
abline(0,0)

dev.off()


jpeg("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/Paper/suppl_figure_2b.jpeg",width=800,height=800,quality=100)

par(mar=c(4,5,3,3))
x <- Rt_wo[d_inc_chd_10$incchd==0]
y <- Rt_w[d_inc_chd_10$incchd==0]-Rt_wo[d_inc_chd_10$incchd==0]


r <- lowess(x[x<0.5],y[x<0.5],f=0.6)
plot(x,y,col="antiquewhite4",xlim=c(0,0.8),ylim=c(-0.10,+0.10),ylab="",xlab="",axes=F,)
title("Non-events",xlab="10-year risk predicted by model containing only FHS risk factors",ylab="Differences in 10-year risk between model \n with FHS risk factors + overall MGRS and only FHS", cex.lab=0.85)
lines(r,lwd=3)
axis(1,at=c(0,0.05,0.10,0.15,0.20,0.30,0.40,0.50,0.60,0.70,0.80),label=c(0,5,10,15,20,30,40,50,60,70,80))
axis(2,at=c(-0.10,-0.05,0,0.05,0.10),label=c(-10,-5,0,+5,+10),las=1)
abline(0,0)

dev.off()


##################################
#### c SUPPLEMENTARY FIGURE 3 ####
##################################

jpeg("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/Paper/suppl_figure_3.jpeg",width=800,height=800,quality=100)

vect <- c(R2_tc,R2_smoke,R2_sbp,R2_bmi,R2_hdl,R2_diab)
names(vect) <- c("TC","SMOKE","SBP","BMI","HDL","DIABETES")
vect_o <- vect[order(vect,decreasing=T)]
barplot(vect_o,space=0.1,col="blue",names.arg=names(vect_o),ylab="Total Variance Explained (%)",ylim=c(0,9))
text(0.67,vect_o[1]+0.4,round(vect_o[1],1), cex=1.05)
text(1.75,vect_o[2]+0.4,round(vect_o[2],1), cex=1.05)
text(2.85,vect_o[3]+0.4,round(vect_o[3],1), cex=1.05)
text(4,vect_o[4]+0.4,round(vect_o[4],1), cex=1.05)
text(5,vect_o[5]+0.4,round(vect_o[5],1), cex=1.05)
text(6,vect_o[6]+0.4,round(vect_o[6],1), cex=1.05)
dev.off()


##########################
##########################
### 7. OTHER ANALYSES ####
##########################
##########################


####################################
### a. C-STATISTIC FROM LOGISTIC ###
####################################

d_inc_chd_p <- d_inc_chd[d_inc_chd$study!="twingene",]

mod1 <- lrm (incchd ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd)
a <- predict(mod1)

mod1_a <- lrm (incchd ~  SEX + ALLCAT_chd_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd)
a_a <- predict(mod1_a)


mod2 <- lrm (incchd ~  SEX  + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd)
b <- predict(mod2)


library(pROC)
roca <- roc(d_inc_chd$incchd ~ a)
roca_a <- roc(d_inc_chd$incchd ~ a_a)
rocb <- roc(d_inc_chd$incchd ~ b)
roc.test(roca,rocb)
roc.test(roca_a,rocb)

## C-statistic for a base model

mod1 <- lrm(incchd ~  ALLCAT_weights_cl_wtccc, data=d_inc_chd)
mod2 <- lrm(incchd ~  ALLCAT_chd_weights_cl_wtccc, data=d_inc_chd)


lrm(incchd ~  age + tc, data=d_inc_chd)


######################
### b. NET BENEFIT ###
######################

source("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Pgm/NetBen.r")

m1 <- cph (Surv(survi, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab + age + study, data=d_inc_chd, surv=T)
r1 <- 1-survest(m1,linear.predictors=m1$linear.predictors, times=10)$surv

m2 <- cph (Surv(survi, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab + age + study + STAR_chd, data=d_inc_chd, surv=T)
r2 <- 1-survest(m2,linear.predictors=m2$linear.predictors, times=10)$surv

dd <- data.frame(d_inc_chd$incchd, d_inc_chd$survi)
colnames(dd) <- c("chd","DURATION")

res1 <- computeNB(patientsToTreat=dd[r1>0.2,],year=10,nTot=nrow(d_inc_chd),N.treat=nrow(dd[r1>0.2,]),theta_values=c(0.7,0.8,0.9))

NB1 <- res1$NB

res2 <- computeNB(patientsToTreat=dd[r2>0.2,],year=10,nTot=nrow(d_inc_chd),N.treat=nrow(dd[r2>0.2,]),theta_values=c(0.7,0.8,0.9))

NB2 <- res2$NB

NB2[2]-NB1[2]


RES <- NULL
## Simulations
for (i in 1:500)
{
	set.seed(123+i)
	d_inc_chd_b <- d_inc_chd[sample(nrow(d_inc_chd), replace=T),]
	m1 <- cph (Surv(survi, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab + age + study, data= d_inc_chd_b, surv=T)
r1 <- 1-survest(m1,linear.predictors=m1$linear.predictors, times=10)$surv

m2 <- cph (Surv(survi, incchd) ~  SEX + tc + hdl + smoke + antihyp + sbp + diab + age + study + STAR_chd, data= d_inc_chd_b, surv=T)
r2 <- 1-survest(m2,linear.predictors=m2$linear.predictors, times=10)$surv

dd <- data.frame(d_inc_chd_b$incchd, d_inc_chd_b$survi)
colnames(dd) <- c("chd","DURATION")

res1 <- computeNB(patientsToTreat=dd[r1>0.2,],year=10,nTot=nrow(d_inc_chd),N.treat=nrow(dd[r1>0.2,]),theta_values=0.8)


res2 <- computeNB(patientsToTreat=dd[r2>0.2,],year=10,nTot=nrow(d_inc_chd),N.treat=nrow(dd[r2>0.2,]),theta_values=0.8)

res <- cbind(res1,res2)
RES <- rbind(RES,res)

}

RES_sd <- sd(RES[,20]-RES[,8])

(NB2[2]-NB1[2])+1.96*RES_sd

(NB2[2]-NB1[2])-1.96*RES_sd


########################
#### SECOND REVIEW #####
########################

library(rmeta)

### Create a forest plot! ###

mod0 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="gender",])

mod1 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="harmony",])

mod2 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="octo",])

mod3 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="satsa",])

mod4 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="twingene",])

mod5 <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="ulsam",])


mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab + study, data=d_inc_chd)


pointest=c(as.numeric(coefficients(mod0)[2]),as.numeric(coefficients(mod1)[2]),as.numeric(coefficients(mod2)[2]),as.numeric(coefficients(mod3)[2]),as.numeric(coefficients(mod4)[2]),as.numeric(coefficients(mod5)[1]))

standest=c(as.numeric(summary(mod0)$coefficients[2,3]),as.numeric(summary(mod1)$coefficients[2,3]),as.numeric(summary(mod2)$coefficients[2,3]),as.numeric(summary(mod3)$coefficients[2,3]),as.numeric(summary(mod4)$coefficients[2,3]),as.numeric(summary(mod5)$coefficients[1,3]))

nnn=c(mod0$nevent,mod1$nevent,mod2$nevent,mod3$nevent,mod4$nevent,mod5$nevent)

# Plot (Supplementary Figure 4)
metaplot(mn=pointest,se=standest,nn=nnn,labels=c("Gender","Hermony","Octo","Satsa","TwinGene","Ulsam"),xlab="Beta for 1 Risk Allele Increase")


### Check association of 100 random score ###

mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_inc_chd[d_inc_chd$study=="twingene",])

d_inc_chd$GWAS_ID <- as.character(d_inc_chd$GWAS_ID)

sc <- read.csv("/Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Data/From_plink/100_scores.csv",stringsAsFactor=F)

d_sc <- merge(d_inc_chd,sc,by.x="GWAS_ID", by.y="V1")

RES <- NULL
for (i in 109:208)
{
	mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + as.numeric(d_sc[,i]) + tc + hdl + smoke + antihyp + sbp + diab, data=d_sc)
	res <- c(as.numeric(coefficients(mod)[2]),as.numeric(summary(mod)$coefficients[2,3]))
	RES <- rbind(RES,res)
}

## Plot (Supplementary Figure 5)

d <- density(RES[,1])
plot(NA, xlab="", ylab="", main="", xlim=c(-0.1,0.3), ylim=c(0,7.5))
abline(v=mean(RES[,1]), lwd=2)
abline(v=0.18691, col="red", lwd=2)
for (i in 1:length(RES[,1]))
{
	abline(v=RES[i,1], col="grey",lwd=0.5)
}
par(new=T)
plot(d, xlab="Beta for 1 Risk Allele Increase", main="",ylab="Density",xlim=c(-0.1,0.3), ylim=c(0,7.5))



### TO calculate the differences we bootstrap


DIFF <- NULL
for (k in 1:200)
{
	set.seed(123+k)
	d_sc_b <- d_sc[sample(nrow(d_sc),replace=T),]
	mod_b <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + ALLCAT_weights_cl_wtccc + tc + hdl + smoke + antihyp + sbp + diab, data=d_sc_b)
	RES <- NULL
	for (i in 109:208)
	{
	mod <- coxph (Surv(age_entry_chd,age_exit_chd, incchd) ~  SEX + as.numeric(d_sc_b[,i]) + tc + hdl + smoke + antihyp + sbp + diab, data=d_sc_b)
	RES <- c(RES,as.numeric(coefficients(mod)[2]))
	}
	DIFF <- c(DIFF,coefficients(mod_b)[2]-mean(RES))
	print(k)
}

1-pnorm(mean(DIFF)/sd(DIFF))

