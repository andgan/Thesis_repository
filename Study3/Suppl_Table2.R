#----------------------------------------------
# Filename: Suppl Table2.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 24SEP2014 - Do vFDR and RDR calculations
# Purpose: Analysis used in the paper. Suppl Table 2
# Note: Main association between each feature in ULSAM and replication in TwinGene, then meta-analysis.
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small)
# Data created: chd_inc_ulsam_PAPER.Rdata chd_inc_twge_PAPER.Rdata chd_inc_meta_PAPER.Rdata Suppl_Table2.csv
#-----------------------------------------------
# OP: R 2.13.1, speedglm
#-----------------------------------------------*/


#### LOAD ULSAM ####

library(survival)

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/ulsam.metabolomics/data_processed/plasma/ulsam_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]


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



# CHD-free subjects
metabo_sub_nopchdT <- metabo_sub_chd[metabo_sub_chd$chd==0 | (metabo_sub_chd$chd==1 & metabo_sub_chd$chdd > metabo_sub_chd$check_d),]

## Log-crp
metabo_sub_nopchdT$logcrp <- log(metabo_sub_nopchdT$crp)
metabo_sub_nopchdT$logtg <- log(metabo_sub_nopchdT$tg)


## Number of metbolites
nnro <- system("awk '{ print NF+1 }' /proj/b2011036/ulsam.metabolomics/data_processed/plasma/ulsam_small.txt", intern=T)
nnro <- as.numeric(nnro)[1]


#### ASSOCIATION ANALYSIS IN ULSAM ###
library(survival)

res <- matrix(NA,nrow=(nnro-2),ncol=4)
for (i in 3:(nnro))
	{
mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU[,i])), data =metabo_sub_nopchdU)
propP <- cox.zph(mod)$table[2,3]

	res[(i-2),] <- c(colnames(metabo_sub_nopchdU)[i],mod$coefficients[2],summary(mod)$coefficients[2,3], propP)		
	print(i)
	}

chd_incS_ulsam10 <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_ulsam10) <- c("feature","beta","se","propTest","p")

## Save results
save(chd_incS_ulsam10,file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/chd_inc_ulsam_PAPER.Rdata")



#### For those with lack of proportionality, rerun it adding a interaction term ####

### Count number of records with lack of proportionality ###
unproptest <- chd_incS_ulsam10[p.adjust(as.numeric(chd_incS_ulsam10[,4]),method="BH")<0.10,]

#### NONE OF THEM LACK OF PROPORTIONALITY #### 


#### REPLICATE FINDINGS (15% FDR in TWINGENE)

chd_incS_ulsam10S <- chd_incS_ulsam10[p.adjust(as.numeric(chd_incS_ulsam10[,5]),method="BH")<0.15,]
### Which P-value correspond ? ###
max(as.numeric(chd_incS_ulsam10[p.adjust(as.numeric(chd_incS_ulsam10[,5]),method="BH")<0.05,5]))


#chd_incS_ulsam10[which(abs(p.adjust(as.numeric(chd_incS_ulsam10[,5]),method="BH")-0.05)==min(abs(p.adjust(as.numeric(chd_incS_ulsam10[,5]),method="BH")-0.05))),]

#### Now I proceed to manual grouping and annotation ####
######## (Andrea_Priority_Table_1_mod_andrea) ########
#### Features to be replicated in Ulsam and TwinGene #####

to_repU <- c("M393.238T397.936","M225.035T32.769","M520.340T357.723","M221.081T145.997","M661.528T590.990","M522.356T395.236","M552.403T462.976","M570.356T377.717","M542.325T331.776","M742.574T674.788","M518.325T328.219","M301.119T163.535","M766.574T661.837","M718.540T633.598","M237.222T457.845","M317.248T445.940","M673.528T584.147","M180.066T113.964","M730.538T628.356","M770.606T700.547","M466.330T470.538","M494.325T336.393","M229.217T385.729","M548.372T413.531","M675.544T608.710","M544.340T351.293","M732.554T648.330","M633.256T229.527","M715.575T635.238","M615.246T242.964","M297.279T516.088","M687.543T602.487","M689.559T624.598")


to_repT <- c("M393.231T384.578","M225.035T32.546","M520.340T346.671","M221.081T142.576","M661.528T580.694","M522.356T382.722","M552.403T449.740","M570.356T366.208","M542.325T321.861","M742.575T667.347","M518.325T318.128","M301.119T159.030","M766.574T654.279","M718.539T624.881","M237.222T443.092","M317.248T431.120","M673.528T573.444","M180.066T112.399","M730.538T619.424","M770.606T693.654","M466.330T457.464","M494.325T325.856","M229.217T373.363","M548.372T400.956","M675.544T599.062","M544.340T340.663","M732.554T640.164","M633.256T223.360","M715.575T626.520","M615.245T236.760","M297.279T501.891","M687.544T592.512","M689.560T615.637")


metabo_sub_nopchdT2 <- metabo_sub_nopchdT[,c(colnames(metabo_sub_nopchdT)[colnames(metabo_sub_nopchdT)%in%to_repT],"age","sex","stratum","incchd","surv_chd","strata1_n","strata2_n","strata3_n","strata4_n")]

metabo_sub_nopchdT2$stratum <- ifelse(metabo_sub_nopchdT2$incchd==1,0,metabo_sub_nopchdT2$stratum)

metabo_sub_nopchdT2$weights <- ifelse(metabo_sub_nopchdT2$stratum==0,1,
ifelse(metabo_sub_nopchdT2$stratum==1,metabo_sub_nopchdT2$strata1_n[1]/sum(metabo_sub_nopchdT2$stratum==1),
ifelse(metabo_sub_nopchdT2$stratum==2,metabo_sub_nopchdT2$strata2_n[1]/sum(metabo_sub_nopchdT2$stratum==2),
ifelse(metabo_sub_nopchdT2$stratum==3,metabo_sub_nopchdT2$strata3_n[1]/sum(metabo_sub_nopchdT2$stratum==3),metabo_sub_nopchdT2$strata4_n[1]/sum(metabo_sub_nopchdT2$stratum==4)))))

### Simple ###
strata <- as.numeric(as.factor(metabo_sub_nopchdT2$weights))-1

res <- NULL
for (i in to_repT)
	{
		stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(metabo_sub_nopchdT2[,i])), data =metabo_sub_nopchdT2,weights=metabo_sub_nopchdT2$weights) 
		dfb<-as.matrix(resid(stratcox,type="dfbeta"))

		# Recalculate the correct SE
		strata_na <- as.numeric(as.factor(metabo_sub_nopchdT2$weights))-1

		gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
		for (s in 1:4)
			{
					indst<-(1:length(metabo_sub_nopchdT2$surv_chd))[strata_na==s]
					m <- as.numeric(table(metabo_sub_nopchdT2$weights)[s+1])
					n <- as.numeric(table(metabo_sub_nopchdT2$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT2$weights)[s+1]))
					if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
			} 
		adjvar<-stratcox$var+gamma
		adjse <- sqrt(diag(adjvar))
		res <- rbind(res,c(i,stratcox$coefficients[3],adjse[3]))			
		print(i)
	}		

	
chd_incS_twge <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_twge) <- c("feature","beta","se","p")

## Save results
save(chd_incS_twge,file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/chd_inc_twge_PAPER.Rdata")


##### META-ANALYZE #####


ulsam_s <- chd_incS_ulsam10S[chd_incS_ulsam10S[,1]%in%to_repU,]
# Order to be in the same order
ulsam_ss <- ulsam_s[order(match(ulsam_s[,1],to_repU)),]
twge_ss <- chd_incS_twge[order(match(chd_incS_twge[,1],to_repU)),]


library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:length(ulsam_s[,1]))
{
	beta <- c(as.numeric(twge_ss[i,2]),as.numeric(ulsam_ss[i,2]))
	se <- c(as.numeric(twge_ss[i,3]),as.numeric(ulsam_ss[i,3]))
	temp <- metagen(beta,se,studlab=c("TwinGene","Ulsam"))
	meta_b <- c(meta_b,temp$TE.random)
	meta_s <- c(meta_s,temp$seTE.random)
}



### Meta-analysis ###
chd_meta <- cbind(ulsam_ss[,1], ulsam_ss[,2],ulsam_ss[,5],twge_ss[,1], twge_ss[,2], twge_ss[,4], meta_b,
2*pnorm(-abs(meta_b/meta_s)))

colnames(chd_meta) <- c("Feature Ulsam","Beta Ulsam","P-value Ulsam","Feature TwinGene","Beta TwinGene","P-value TwinGene","beta meta","p meta")

write.csv(chd_meta, file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Suppl_Table2.csv", quote=F, row.names=F)


### Check sign test ###
binom.test(sum(sign(as.numeric(chd_meta[,2]))==sign(as.numeric(chd_meta[,5]))), nrow(chd_meta), alternative="greater") 


binom.test(27, 33, alternative="greater") 



####################################################################
###### REVIEWER ANALYSIS: ANALYSIS TO ESTIMATED RDR AND vFDR #######
####################################################################


load("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/chd_inc_ulsam_PAPER.Rdata")


ct<-	qt(10^(-3.10)/2,df=131+897-2,lower.tail=FALSE)	
cv<-	qt(10^(-1.3)/2,df=282+1388-2,lower.tail=FALSE)

RDR.est<-rdr.est(as.numeric(chd_incS_ulsam10[,2])/as.numeric(chd_incS_ulsam10[,3]), nt0=131, nt1=897, nv0=282, nv1=1388,  c.t=3.364808, c.v=1.960371)

qvalue(as.numeric(chd_incS_ulsam10[,5]))$qvalues[chd_incS_ulsam10[,1]%in%to_repU]



###############################################################################
###### NOT INCLUDED IN THE PAPER - ANALYSIS USING CLUSTERING ALGHORITHM #######
###############################################################################

## Load data ##
load("/lynx/cvol/v38/b2011036/metabolomics/ulsam_plasma_ramclust_nosing.RData")
load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")
load("/lynx/cvol/v38/b2011036/metabolomics/ulsam_plasma_ramclust_withsing.RData")


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


## KEEP ONLY USEFUL VARIABLES ##
metabo_sub_nopchdUB <- metabo_sub_nopchdU[,c("pat_time","pat","diab","tc","hdl","age","sbp","smoke","bmi","incchd10","surv_chd10")]
clu <- data.frame(ulsam_plasma_ramclust_no_singletons$SpecAbundAve, stringAsFactors=F)
clu$id <- rownames(clu)

## MERGE ##
cluUB <- merge(clu, metabo_sub_nopchdUB, by.x="id", by.y="pat_time")


nnroCLU <- ncol(ulsam_plasma_ramclust_no_singletons$SpecAbundAve)
library(survival)

res <- matrix(NA,nrow=(nnroCLU-2),ncol=4)
for (i in 2:(nnroCLU-1))
	{
mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(cluUB[,i])), data=cluUB)
propP <- cox.zph(mod)$table[2,3]

	res[(i-1),] <- c(colnames(cluUB)[i],mod$coefficients[2],summary(mod)$coefficients[2,3], propP)		
	print(i)
	}

chd_incS_ulsam10CLU <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_ulsam10CLU) <- c("feature","beta","se","propTest","p")

## DO THE SAME INCLUDING ALSO SINGLE CLUSTERS ##
clu2 <- data.frame(ulsam_plasma_ramclust_with_singletons$SpecAbundAve, stringAsFactors=F)
clu2$id <- rownames(clu2)


cluUB2 <- merge(clu2, metabo_sub_nopchdUB, by.x="id", by.y="pat_time")


nnroCLU2 <- ncol(ulsam_plasma_ramclust_with_singletons$SpecAbundAve)
library(survival)

res <- matrix(NA,nrow=(nnroCLU2-2),ncol=4)
for (i in 2:(nnroCLU2-1))
	{
mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(cluUB2[,i])), data=cluUB2)
propP <- cox.zph(mod)$table[2,3]

	res[(i-1),] <- c(colnames(cluUB2)[i],mod$coefficients[2],summary(mod)$coefficients[2,3], propP)		
	print(i)
	}

chd_incS_ulsam10CLU2 <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_ulsam10CLU2) <- c("feature","beta","se","propTest","p")


### NOW COMPARE THE RESULTS ###
single <- chd_incS_ulsam10[p.adjust(as.numeric(chd_incS_ulsam10[,5]),method="BH")<0.02,]
group <- chd_incS_ulsam10CLU[p.adjust(as.numeric(chd_incS_ulsam10CLU[,5]),method="BH")<0.02,]
group2 <- chd_incS_ulsam10CLU2[p.adjust(as.numeric(chd_incS_ulsam10CLU2[,5]),method="BH")<0.02,]


write.csv(single,file="test.csv")
write.csv(group,file="test2.csv")
write.csv(group2,file="test3.csv")



#################################################################
###### NOT INCLUDED IN THE PAPER - ANALYSIS USING 10% FDR #######
#################################################################


library(survival)

res <- matrix(NA,nrow=length(com_PTU[!is.na(com_PTU[,2]) & !is.na(com_PTU[,3]),3]),ncol=4)
for (i in 1:length(com_PTU[!is.na(com_PTU[,2]) & !is.na(com_PTU[,3]),3]))
	{
mod <- coxph(Surv(surv_chd10,incchd10) ~ age+scale(as.numeric(metabo_sub_nopchdU2[,i])), data =metabo_sub_nopchdU2)
propP <- cox.zph(mod)$table[2,3]

	res[i,] <- c(colnames(metabo_sub_nopchdU2)[i],mod$coefficients[2],summary(mod)$coefficients[2,3], propP)		
	print(i)
	}

chd_incS_ulsam102 <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_ulsam102) <- c("feature","beta","se","propTest","p")


### If you want to use 10%FDR ###
metabo_sub_nopchdT2 <- metabo_sub_nopchdT[,c(colnames(metabo_sub_nopchdT)[colnames(metabo_sub_nopchdT)%in%com_PTU[!is.na(com_PTU[,2]) & !is.na(com_PTU[,3]),2]],"age","sex","stratum","incchd","surv_chd","strata1_n","strata2_n","strata3_n","strata4_n")]
commonS <- com_PTU[!is.na(com_PTU[,2]) & !is.na(com_PTU[,3]),]



metabo_sub_nopchdT2 <- metabo_sub_nopchdT

metabo_sub_nopchdT2$stratum <- ifelse(metabo_sub_nopchdT2$incchd==1,0,metabo_sub_nopchdT2$stratum)

metabo_sub_nopchdT2$weights <- ifelse(metabo_sub_nopchdT2$stratum==0,1,
ifelse(metabo_sub_nopchdT2$stratum==1,metabo_sub_nopchdT2$strata1_n[1]/sum(metabo_sub_nopchdT2$stratum==1),
ifelse(metabo_sub_nopchdT2$stratum==2,metabo_sub_nopchdT2$strata2_n[1]/sum(metabo_sub_nopchdT2$stratum==2),
ifelse(metabo_sub_nopchdT2$stratum==3,metabo_sub_nopchdT2$strata3_n[1]/sum(metabo_sub_nopchdT2$stratum==3),metabo_sub_nopchdT2$strata4_n[1]/sum(metabo_sub_nopchdT2$stratum==4)))))

### Simple ###
strata <- as.numeric(as.factor(metabo_sub_nopchdT2$weights))-1

res <- matrix(NA,nrow=9754-3,ncol=3)
for (i in 3:9754)
	{
		stratcox<-coxph(Surv(surv_chd,incchd) ~ age + sex + scale(as.numeric(metabo_sub_nopchdT2[,i])), data =metabo_sub_nopchdT2,weights=metabo_sub_nopchdT2$weights) 
		dfb<-as.matrix(resid(stratcox,type="dfbeta"))

		# Recalculate the correct SE
		strata_na <- as.numeric(as.factor(metabo_sub_nopchdT2$weights))-1

		gamma<-matrix(0,dim(stratcox$var)[1],dim(stratcox$var)[2])
		for (s in 1:4)
			{
					indst<-(1:length(metabo_sub_nopchdT2$surv_chd))[strata_na==s]
					m <- as.numeric(table(metabo_sub_nopchdT2$weights)[s+1])
					n <- as.numeric(table(metabo_sub_nopchdT2$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT2$weights)[s+1]))
					if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
			} 
		adjvar<-stratcox$var+gamma
		adjse <- sqrt(diag(adjvar))
		res[i-3,] <- c(colnames(metabo_sub_nopchdT2)[i],stratcox$coefficients[3],adjse[3])			
		print(i)
	}		

	
chd_incS_twge <- cbind(res,2*pnorm(-abs(as.numeric(res[,2])/as.numeric(res[,3]))))
colnames(chd_incS_twge) <- c("feature","beta","se","p")

chd_incS_twge[chd_incS_twge[,1]=="M746.569T654.597"]



ind <- which(!chd_incS_ulsam102[,1]%in%ulsam_ss[,1])


### Meta-analysis ###
chd_meta <- cbind(c(ulsam_ss[,1], chd_incS_ulsam102[ind,1]),
c(ulsam_ss[,2],chd_incS_ulsam102[ind,2]),
c(ulsam_ss[,5],chd_incS_ulsam102[ind,5]),
c(twge_ss[,1],rep("-",length(ind))),
c(twge_ss[,2],rep("-",length(ind))),
c(twge_ss[,4],rep("-",length(ind))),
c(meta_b,rep("-",length(ind))),
c(2*pnorm(-abs(meta_b/meta_s)),rep("-",length(ind))))

colnames(chd_meta) <- c("Feature Ulsam","Beta Ulsam","P-value Ulsam","Feature TwinGene","Beta TwinGene","P-value TwinGene","beta meta","p meta")

## Features with a P-value < 0.05 in TwinGene
chd_meta <- cbind(chd_meta,ifelse(as.numeric(chd_meta[,6])<0.05,1,0))
colnames(chd_meta)[9] <- c("Sign in TwinGene")

#save(chd_meta,file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/chd_inc_meta_PAPER.Rdata")
write.csv(chd_meta, file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/material_for_Suppl_Table2_test.csv", quote=F, row.names=F)

chd_metaM <- chd_meta[chd_meta[,8]!="-",]
unproptest <- chd_metaM[p.adjust(as.numeric(chd_m
	


