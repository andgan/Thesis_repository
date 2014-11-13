#----------------------------------------------
# Filename: Figure2_PanelC_suppl_fig_4.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 27FEB2014
# Updated: 12JUL2014 - Update figures names and covariate adjustments
# Purpose: Analysis used in the paper. Figure2_PanelC and SUPPLEMENTARY FIGURE  4
# Note: Association between CHD SNPs and 4 main metabolites. Note that it uses SNPs extracted by extract_snps.py
#  !!! In the end of the script there is also the procedure to check for rare variants in some candidate genes using the ecxome chip, it does not mean that we will include this analysis in the final paper, since the results are weak.
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small, Pivus small)  exp_metabo_pivus.csv exp_metabo_ulsam.csv exp_metabo_twge.csv Annotated_snplist.txt
# Data created: Fig_2_p3.pdf Suppl_Figure4.pdf
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/



############
## PIVUS ##
###########

load("/home/andrea/glob/alignment_pivus_small/Results/Final_datasets/step4.Rdata")

## 
# M494.325T326.430 = LysoPC 16:1
# M496.340T365.792 = LysoPC 16:0
# M522.356T384.872 = LysoPC 18:1
# M524.372T425.563 = LysoPC 18:0
# M520.340T347.575 = LysoPC 18:2
# M566.322T352.005 = LysoPC 20:4

# M730.538T619.422 = PC(32:2)
# M734.560T640.059 = PC(32:0)
# M812.616T686.543 = PC(36:0)
# M319.195T417.146 = LIONELIC ACID 


# M760.585T667.435 = PC(34:1)
# M758.569T649.183 = PC(34:2)
# M788.616T693.450 = PC(36:1)
# M786.601T676.687== PC(36:2)


# M118.087T33.174 = Betaine
# M105.111T32.218 = Choline

# M377.267T388.196 = MG18:2
# M661.528T579.938 = PE_cer

mPIVUS <- metabo_p[,c("id","M494.325T326.430","M496.340T365.792","M522.356T384.872","M524.372T425.563","M520.340T347.575","M566.322T352.005",
"M730.538T619.422","M734.560T640.059","M786.601T676.687","M812.616T686.543",
"M319.195T417.146",
"M118.087T33.174","M105.111T32.218",
"M760.585T667.435","M758.569T649.183","M788.616T693.450","M377.267T388.196","M661.528T579.938","check_d","check_d","smoke01","diab","tg","antihyp")]
cardiphen <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS pek tander vs athero6.txt", header=T, sep="\t", stringsAsFactor=F)
pla2 <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS data PLA2s.txt", header=T, sep="\t", stringsAsFactor=F)



# Numeric
mPIVUS$id <- as.numeric(mPIVUS$id)
pla2 <- apply(pla2,2,as.numeric)

# Merging
tmp <- merge(mPIVUS,cardiphen,by="id")
mPIVUSF <- merge(tmp,pla2,by.x="id", by.y="lpnr")

mPIVUSF$logcrp <- log(mPIVUSF$crp)
mPIVUSF$logtg <- log(mPIVUSF$tg)


colnames(mPIVUSF)[2:19] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC32_2","PC32_0","PC36_2","PC36_0",
"LIONELIC",
"BETAINE","CHOLINE",
"PC34_1","PC34_2","PC36_1","MG18_2","PE_cer")


mPIVUSF$gwas_id <- paste("PIVUS",sprintf( "%04d", as.numeric(mPIVUSF$id) ),sep="")

### Load Extracted data ###
psnp <- read.csv("/home/andrea/glob/alignment_pivus_small/Data/exp_metabo_pivus.csv")

### Merge with data ####
mPIVUSF2 <- merge(psnp,mPIVUSF,by.x="ID", by.y="gwas_id")

#temp181a <- (mPIVUSF2$LPC18_1~mPIVUSF2[,i] + logcrp+manuelltsbp+ hdl + bmi + smoke01 + diab+ kolesterol+kvinna1,data=mPIVUSF2)

PRES181 <- NULL
PRES182 <- NULL
PRESmg <- NULL
PRESPE_cer <- NULL
PRESmga <- NULL
PRESPE_cera <- NULL
PRES181a <- NULL
PRES182a <- NULL

for(i in colnames(psnp)[1:(ncol(psnp)-1)])
{
	temp181 <- lm(scale(as.numeric(mPIVUSF2$LPC18_1))~mPIVUSF2[,i] + kvinna1,data=mPIVUSF2)
	temp182 <- lm(scale(as.numeric(mPIVUSF2$LPC18_2))~mPIVUSF2[,i] + kvinna1,data=mPIVUSF2)
	tempmg <-  lm(scale(as.numeric(mPIVUSF2$MG18_2))~mPIVUSF2[,i] + kvinna1,data=mPIVUSF2)
	tempPE_cer <- lm(scale(as.numeric(mPIVUSF2$PE_cer))~mPIVUSF2[,i] + kvinna1,data=mPIVUSF2)
	tempmga <- lm(scale(as.numeric(mPIVUSF2$MG18_2))~mPIVUSF2[,i] +manuelltsbp+ hdl + bmi + smoke01 + diab+ ldl + antihyp+logtg,data=mPIVUSF2)
	tempPE_cera <- lm(scale(as.numeric(mPIVUSF2$PE_cer))~mPIVUSF2[,i] +manuelltsbp+ hdl + bmi + smoke01 + diab+ ldl + antihyp+logtg,data=mPIVUSF2)
	temp181a <- lm(scale(as.numeric(mPIVUSF2$LPC18_1))~mPIVUSF2[,i] + kvinna1 +manuelltsbp+ hdl + bmi + smoke01 + diab+ ldl + antihyp+logtg,data=mPIVUSF2)
	temp182a <- lm(scale(as.numeric(mPIVUSF2$LPC18_2))~mPIVUSF2[,i] + kvinna1 +manuelltsbp+ hdl + bmi + smoke01 + diab+ ldl + antihyp+logtg,data=mPIVUSF2)
	
	
	PRES181 <- rbind(PRES181,c(i,temp181$coefficients[2],summary(temp181)$coefficients[2,2]))
	PRES182 <- rbind(PRES182,c(i,temp182$coefficients[2],summary(temp182)$coefficients[2,2]))
	PRESmg <- rbind(PRESmg,c(i,tempmg$coefficients[2],summary(tempmg)$coefficients[2,2]))
	PRESPE_cer <- rbind(PRESPE_cer,c(i,tempPE_cer$coefficients[2],summary(tempPE_cer)$coefficients[2,2]))
	PRESmga <- rbind(PRESmga,c(i,tempmga$coefficients[2],summary(tempmga)$coefficients[2,2]))
	PRESPE_cera <- rbind(PRESPE_cera,c(i,tempPE_cera$coefficients[2],summary(tempPE_cera)$coefficients[2,2]))
	PRES181a <- rbind(PRES181a,c(i,temp181a$coefficients[2],summary(temp181a)$coefficients[2,2]))
	PRES182a <- rbind(PRES182a,c(i,temp182a$coefficients[2],summary(temp182a)$coefficients[2,2]))
}


colnames(PRES181) <- c("snp","Pb","Pse")
colnames(PRES182) <- c("snp","Pb","Pse")
colnames(PRESmg) <- c("snp","Pb","Pse")
colnames(PRESPE_cer) <- c("snp","Pb","Pse")
colnames(PRESmga) <- c("snp","Pb","Pse")
colnames(PRESPE_cera) <- c("snp","Pb","Pse")
colnames(PRES181a) <- c("snp","Pb","Pse")
colnames(PRES182a) <- c("snp","Pb","Pse")


##############
## TWINGENE ##
##############

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

# M393.231T384.578 = MG18:2
# M661.528T580.694 = PE_cer

mTWINGENEF <- metabo_p[,c("twinnr","M494.325T325.856","M496.340T364.064","M522.356T382.722","M524.372T423.532","M520.340T346.671","M566.322T340.659",
"M730.538T619.424","M732.554T640.164","M734.560T640.164","M782.570T652.163",
"M786.601T677.527","M812.616T687.249",
"M285.279T515.220","M283.264T452.573","M321.270T463.343","M551.434T435.908","M293.179T392.421","M281.010T31.737","M319.195T415.326",
"M118.087T33.403","M105.111T32.168",
"M760.585T667.863","M758.570T649.594","M788.617T694.035","M393.231T384.578","M661.528T580.694",colnames(metabo_p)[9756:ncol(metabo_p)])]

colnames(mTWINGENEF)[2:27] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC32_2","PC32_1","PC32_0","PC36_4","PC36_2","PC38_3",
"STEARIC","OLEIC","VACCENIC","PALMITIC","PALMITOLEIC","GLYPC","LIONELIC",
"BETAINE","CHOLINE",
"PC34_1","PC34_2","PC36_1","MG18_2","PE_cer")

mTWINGENEF$logcrp <- log(mTWINGENEF$crp)
mTWINGENEF$logtg <- log(mTWINGENEF$tg)
### Recode diabetes ###
mTWINGENEF$diab_baseline_new <- ifelse((mTWINGENEF$diabetess==1 & mTWINGENEF$diabetess_doctor==1) | mTWINGENEF$diabetess_medicine==1 | (mTWINGENEF$glucos>= 7 & mTWINGENEF$fasting==0) | (mTWINGENEF$glucos>= 11 & mTWINGENEF$fasting!=0) | (mTWINGENEF$diab_regd <= mTWINGENEF$check_d & !is.na(mTWINGENEF$diab_regd)),1,0)


### Load Extracted data ###
psnt <- read.csv("/home/andrea/glob/alignment_twge_small/Data/exp_metabo_twge.csv", stringsAsFactor=F)

### Merge with data ####
mTWINGENEF2 <- merge(psnt,mTWINGENEF,by.x="ID", by.y="gwas_id")


TRES181 <- NULL
TRES182 <- NULL
TRESmg <- NULL
TRESPE_cer <- NULL
TRESmga <- NULL
TRESPE_cera <- NULL
TRES181a <- NULL
TRES182a <- NULL

for(i in colnames(psnt)[1:(ncol(psnt)-1)])
{
	temp181 <- lm(scale(as.numeric(mTWINGENEF2$LPC18_1))~mTWINGENEF2[,i] + age + sex,data=mTWINGENEF2)
	temp182 <- lm(scale(as.numeric(mTWINGENEF2$LPC18_2))~mTWINGENEF2[,i] + age + sex,data=mTWINGENEF2)
	tempmg <-  lm(scale(as.numeric(mTWINGENEF2$MG18_2))~mTWINGENEF2[,i] + age + sex,data=mTWINGENEF2)
	tempPE_cer <- lm(scale(as.numeric(mTWINGENEF2$PE_cer))~mTWINGENEF2[,i] + age + sex,data=mTWINGENEF2)
	tempmga <- lm(scale(as.numeric(mTWINGENEF2$MG18_2))~mTWINGENEF2[,i] + age+sex+sbp+ hdl + bmi + smoke01 + diab_baseline+ ldl+antihyp + logtg,data=mTWINGENEF2)
	tempPE_cera <- lm(scale(as.numeric(mTWINGENEF2$PE_cer))~mTWINGENEF2[,i] + age+sex+sbp+ hdl + bmi + smoke01 + diab_baseline+ tc+antihyp + logtg,data=mTWINGENEF2)
	temp181a <- lm(scale(as.numeric(mTWINGENEF2$LPC18_1))~mTWINGENEF2[,i] + age+sex+sbp+ hdl + bmi + smoke01 + diab_baseline+ ldl+antihyp + logtg,data=mTWINGENEF2)
	temp182a <- lm(scale(as.numeric(mTWINGENEF2$LPC18_2))~mTWINGENEF2[,i] + age+sex+sbp+ hdl + bmi + smoke01 + diab_baseline+ ldl+antihyp + logtg,data=mTWINGENEF2)
	
	
	TRES181 <- rbind(TRES181,c(i,temp181$coefficients[2],summary(temp181)$coefficients[2,2]))
	TRES182 <- rbind(TRES182,c(i,temp182$coefficients[2],summary(temp182)$coefficients[2,2]))
	TRESmg <- rbind(TRESmg,c(i,tempmg$coefficients[2],summary(tempmg)$coefficients[2,2]))
	TRESPE_cer <- rbind(TRESPE_cer,c(i,tempPE_cer$coefficients[2],summary(tempPE_cer)$coefficients[2,2]))
	TRESmga <- rbind(TRESmga,c(i,tempmga$coefficients[2],summary(tempmga)$coefficients[2,2]))
	TRESPE_cera <- rbind(TRESPE_cera,c(i,tempPE_cera$coefficients[2],summary(tempPE_cera)$coefficients[2,2]))
	TRES181a <- rbind(TRES181a,c(i,temp181a$coefficients[2],summary(temp181a)$coefficients[2,2]))
	TRES182a <- rbind(TRES182a,c(i,temp182a$coefficients[2],summary(temp182a)$coefficients[2,2]))
}


colnames(TRES181) <- c("snp","Tb","Tse")
colnames(TRES182) <- c("snp","Tb","Tse")
colnames(TRESmg) <- c("snp","Tb","Tse")
colnames(TRESPE_cer) <- c("snp","Tb","Tse")
colnames(TRESmga) <- c("snp","Pb","Pse")
colnames(TRESPE_cera) <- c("snp","Tb","Tse")
colnames(TRES181a) <- c("snp","Tb","Tse")
colnames(TRES182a) <- c("snp","Tb","Tse")


###########
## ULSAM ##
###########


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

# M393.238T397.936= MG18:2
# M661.528T590.990 = PE_cer

mULSAMF <- metabo_p[metabo_p$time==0 & metabo_p$double_==0,c("pat_time","M494.324T324.997","M496.340T375.465","M522.356T395.236","M524.372T436.589","M520.340T357.723","M566.322T351.269",
"M730.538T628.356","M732.554T648.330","M734.560T648.335","M782.575T647.537","M786.609T683.883","M812.616T693.845",
"M285.280T528.734","M283.264T466.637","M321.270T477.975","M551.434T449.611","M293.179T405.566","M281.010T31.036","M319.195T428.855",
"M118.087T34.235","M105.111T32.235",
"M760.594T674.770","M758.569T657.575","M788.616T700.812","M393.238T397.936","M661.528T590.990",colnames(metabo_p)[10163:ncol(metabo_p)])]

colnames(mULSAMF)[2:27] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC32_2","PC32_1","PC32_0","PC36_4","PC36_2","PC38_3",
"STEARIC","OLEIC","VACCENIC","PALMITIC","PALMITOLEIC","GLYPC","LIONELIC",
"BETAINE","CHOLINE",
"PC34_1","PC34_2","PC36_1","MG18_2","PE_cer")


mULSAMF$logcrp <- log(mULSAMF$crp)
mULSAMF$logtg <- log(mULSAMF$tg)
mULSAMF$gwas_id <- paste("ULSAM",sprintf( "%04d", as.numeric(mULSAMF$pat) ),sep="")

### Load Extracted data ###
psnu <- read.csv("/home/andrea/glob/alignment_ulsam_small/Data/exp_metabo_ulsam.csv")

### Merge with data ####
mULSAMF2 <- merge(psnu,mULSAMF,by.x="ID", by.y="gwas_id")


URES181 <- NULL
URES182 <- NULL
URESmg <- NULL
URESPE_cer <- NULL
URESmga <- NULL
URESPE_cera <- NULL
URES181a <- NULL
URES182a <- NULL

for(i in colnames(psnu)[1:(ncol(psnu)-1)])
{
	temp181 <- lm(scale(as.numeric(mULSAMF2$LPC18_1))~mULSAMF2[,i],data=mULSAMF2)
	temp182 <- lm(scale(as.numeric(mULSAMF2$LPC18_2))~mULSAMF2[,i],data=mULSAMF2)
	tempmg <-  lm(scale(as.numeric(mULSAMF2$MG18_2))~mULSAMF2[,i],data=mULSAMF2)
	tempPE_cer <- lm(scale(as.numeric(mULSAMF2$PE_cer))~mULSAMF2[,i],data=mULSAMF2)
	tempmga <- lm(scale(as.numeric(mULSAMF2$MG18_2))~mULSAMF2[,i] +sbp+ hdl + bmi + smoke + diab+ ldl+antihyp+logtg,data=mULSAMF2)
	tempPE_cera <- lm(scale(as.numeric(mULSAMF2$PE_cer))~mULSAMF2[,i]+sbp+ hdl + bmi + smoke + diab+ ldl+antihyp+logtg,data=mULSAMF2)
	temp181a <- lm(scale(as.numeric(mULSAMF2$LPC18_1))~mULSAMF2[,i]+sbp+ hdl + bmi + smoke + diab+ ldl+antihyp+logtg,data=mULSAMF2)
	temp182a <- lm(scale(as.numeric(mULSAMF2$LPC18_2))~mULSAMF2[,i]+sbp+ hdl + bmi + smoke + diab+ ldl+antihyp+logtg,data=mULSAMF2)
	
	
	URES181 <- rbind(URES181,c(i,temp181$coefficients[2],summary(temp181)$coefficients[2,2]))
	URES182 <- rbind(URES182,c(i,temp182$coefficients[2],summary(temp182)$coefficients[2,2]))
	URESmg <- rbind(URESmg,c(i,tempmg$coefficients[2],summary(tempmg)$coefficients[2,2]))
	URESPE_cer <- rbind(URESPE_cer,c(i,tempPE_cer$coefficients[2],summary(tempPE_cer)$coefficients[2,2]))
	URESmga <- rbind(URESmga,c(i,tempmga$coefficients[2],summary(tempmga)$coefficients[2,2]))
	URESPE_cera <- rbind(URESPE_cera,c(i,tempPE_cera$coefficients[2],summary(tempPE_cera)$coefficients[2,2]))
	URES181a <- rbind(URES181a,c(i,temp181a$coefficients[2],summary(temp181a)$coefficients[2,2]))
	URES182a <- rbind(URES182a,c(i,temp182a$coefficients[2],summary(temp182a)$coefficients[2,2]))
}

colnames(URES181) <- c("snp","Ub","Use")
colnames(URES182) <- c("snp","Ub","Use")
colnames(URESmg) <- c("snp","Ub","Use")
colnames(URESPE_cer) <- c("snp","Ub","Use")
colnames(URESmga) <- c("snp","Ub","Use")
colnames(URESPE_cera) <- c("snp","Ub","Use")
colnames(URES181a) <- c("snp","Ub","Use")
colnames(URES182a) <- c("snp","Ub","Use")


#########################
##### META-ANALYSIS #####
#########################


#### RES181 ####

RES181 <- merge(URES181,TRES181, by="snp",all.y=T)
RES181 <- merge(PRES181,RES181, by="snp",all.y=T)


library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RES181))
{
	beta <- c(as.numeric(as.character(RES181[i,2])),as.numeric(as.character(RES181[i,4])),as.numeric(as.character(RES181[i,6])))
	se <- c(as.numeric(as.character(RES181[i,3])),as.numeric(as.character(RES181[i,5])),as.numeric(as.character(RES181[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


meta181 <- cbind(as.character(RES181[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))



#### RES182 ####


RES182 <- merge(URES182,TRES182, by="snp",all.y=T)
RES182 <- merge(PRES182,RES182, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RES182))
{
	beta <- c(as.numeric(as.character(RES182[i,2])),as.numeric(as.character(RES182[i,4])),as.numeric(as.character(RES182[i,6])))
	se <- c(as.numeric(as.character(RES182[i,3])),as.numeric(as.character(RES182[i,5])),as.numeric(as.character(RES182[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


meta182 <- cbind(as.character(RES182[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))



#### RESMG ####

RESmg <- merge(URESmg,TRESmg, by="snp",all.y=T)
RESmg <- merge(PRESmg,RESmg, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RESmg))
{
	beta <- c(as.numeric(as.character(RESmg[i,2])),as.numeric(as.character(RESmg[i,4])),as.numeric(as.character(RESmg[i,6])))
	se <- c(as.numeric(as.character(RESmg[i,3])),as.numeric(as.character(RESmg[i,5])),as.numeric(as.character(RESmg[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


metamg <- cbind(as.character(RESmg[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))




#### PE_cer ####


RESPE_cer <- merge(URESPE_cer,TRESPE_cer, by="snp",all.y=T)
RESPE_cer <- merge(PRESPE_cer,RESPE_cer, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RESPE_cer))
{
	beta <- c(as.numeric(as.character(RESPE_cer[i,2])),as.numeric(as.character(RESPE_cer[i,4])),as.numeric(as.character(RESPE_cer[i,6])))
	se <- c(as.numeric(as.character(RESPE_cer[i,3])),as.numeric(as.character(RESPE_cer[i,5])),as.numeric(as.character(RESPE_cer[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


metaPE_cer <- cbind(as.character(RESPE_cer[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))




#### RESMGA ####

RESmga <- merge(URESmga,TRESmga, by="snp",all.y=T)
RESmga <- merge(PRESmga,RESmga, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RESmg))
{
	beta <- c(as.numeric(as.character(RESmga[i,2])),as.numeric(as.character(RESmga[i,4])),as.numeric(as.character(RESmga[i,6])))
	se <- c(as.numeric(as.character(RESmga[i,3])),as.numeric(as.character(RESmga[i,5])),as.numeric(as.character(RESmga[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


metamga <- cbind(as.character(RESmga[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))




#### RES182/HDLa ####


RESPE_cera <- merge(URESPE_cera,TRESPE_cera, by="snp",all.y=T)
RESPE_cera <- merge(PRESPE_cera,RESPE_cera, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RESPE_cera))
{
	beta <- c(as.numeric(as.character(RESPE_cera[i,2])),as.numeric(as.character(RESPE_cera[i,4])),as.numeric(as.character(RESPE_cera[i,6])))
	se <- c(as.numeric(as.character(RESPE_cera[i,3])),as.numeric(as.character(RESPE_cera[i,5])),as.numeric(as.character(RESPE_cera[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


metaPE_cera <- cbind(as.character(RESPE_cera[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))



#### RES181A ####

RES181a <- merge(URES181a,TRES181a, by="snp",all.y=T)
RES181a <- merge(PRES181a,RES181a, by="snp",all.y=T)


library(meta)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RES181a))
{
	beta <- c(as.numeric(as.character(RES181a[i,2])),as.numeric(as.character(RES181a[i,4])),as.numeric(as.character(RES181a[i,6])))
	se <- c(as.numeric(as.character(RES181a[i,3])),as.numeric(as.character(RES181a[i,5])),as.numeric(as.character(RES181a[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


meta181a <- cbind(as.character(RES181a[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))



#### RES182 ####


RES182a <- merge(URES182a,TRES182a, by="snp",all.y=T)
RES182a <- merge(PRES182a,RES182a, by="snp",all.y=T)

meta_b <- NULL
meta_s <- NULL
for (i in 1:nrow(RES182a))
{
	beta <- c(as.numeric(as.character(RES182a[i,2])),as.numeric(as.character(RES182a[i,4])),as.numeric(as.character(RES182a[i,6])))
	se <- c(as.numeric(as.character(RES182a[i,3])),as.numeric(as.character(RES182a[i,5])),as.numeric(as.character(RES182a[i,7])))
	temp <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))
	meta_b <- c(meta_b,temp$TE.fixed)
	meta_s <- c(meta_s,temp$seTE.fixed)
}


meta182a <- cbind(as.character(RES182a[,1]), meta_b,2*pnorm(-abs(meta_b/meta_s)))






##################################
##### PLOT FIGURE 2- PANEL C #####
##################################

## Read annotation file ##
ann <- read.table("/home/andrea/glob/alignment_pivus_small/Data/Annotated_snplist.txt", stringsAsFactor=F)
snp_meta <- sapply(strsplit(meta181[,1],"_",fixed=T),"[[",1)
FF <- cbind(snp_meta,-log10(as.numeric(meta181[,3])),-log10(as.numeric(meta182[,3])),-log10(as.numeric(metamg[,3])),-log10(as.numeric(metaPE_cer[,3])))
FF <- FF[FF[,1]%in%ann[,2],]
FFo <- FF[order(match(FF[,1],ann[,2])),]


#meta181b <- meta181[snp_meta%in%ann[ann[,5]=="CHD",2],]
#metamgb <- metamg[snp_meta%in%ann[ann[,5]=="CHD",2],]


y <- FFo[,2:ncol(FFo)]
#rownames(y) <- ann[order(match(ann[ann[,5]=="CHD",2],snp_meta[snp_meta%in%ann[ann[,5]=="CHD",2]])),4]
rownames(y) <- ann[,4]
colnames(y) <- c("LysoPC 18:1","LysoPC 18:2","MG 18:2","SM 28:1")

write.csv(y,file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/fig4_noadj.csv")

grid  <- expand.grid(x=seq(1:ncol(y)), y=seq(1:nrow(y)))
grid$z <- as.numeric(t(y))
grid$sign <- ifelse(as.numeric(t(y))>log10((length(y)/0.05)),"","")

library(lattice)
library(gplots)
library(gplots)

pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p3.pdf", height=8, width=5)
levelplot(z~x*y,grid, border="blue",cex.main=0.5, col.regions=rev(colorRampPalette(c("red","green", "aliceblue"), bias=0.2,interpolate="linear")(200)),xlab="",ylab="",at=c(0,1.30,seq(2,9,0.7)), colorkey=list(space="right"),scales=list(x=list(at=seq(1:ncol(y)),labels=colnames(y),rot=45),y = list(at=seq(1:nrow(y)),labels=rownames(y),font=3)),
panel=function(...) {
                       arg <- list(...)
                       panel.levelplot(...,)
                       panel.text(grid$x, grid$y,grid$sign)})
dev.off()



################################
#### SUPPLEMENTARY FIGURE  4 ###
################################



## Read annotation file ##
ann <- read.table("/home/andrea/glob/alignment_pivus_small/Data/Annotated_snplist.txt", stringsAsFactor=F)
snp_meta <- sapply(strsplit(meta181a[,1],"_",fixed=T),"[[",1)
FF <- cbind(snp_meta,-log10(as.numeric(meta181a[,3])),-log10(as.numeric(meta182a[,3])),-log10(as.numeric(metamga[,3])),-log10(as.numeric(metaPE_cera[,3])))
FF <- FF[FF[,1]%in%ann[,2],]
FFo <- FF[order(match(FF[,1],ann[,2])),]


#meta181b <- meta181[snp_meta%in%ann[ann[,5]=="CHD",2],]
#metamgb <- metamg[snp_meta%in%ann[ann[,5]=="CHD",2],]


y <- FFo[,2:ncol(FFo)]
#rownames(y) <- ann[order(match(ann[ann[,5]=="CHD",2],snp_meta[snp_meta%in%ann[ann[,5]=="CHD",2]])),4]
rownames(y) <- ann[,4]
colnames(y) <- c("LysoPC 18:1","LysoPC 18:2","MG 18:2","SM 28:1")

write.csv(y,file="/home/andrea/glob/alignment_pivus_small/Results/chd/paper/fig4_adj.csv")


grid  <- expand.grid(x=seq(1:ncol(y)), y=seq(1:nrow(y)))
grid$z <- as.numeric(t(y))
grid$sign <- ifelse(as.numeric(t(y))>log10((length(y)/0.05)),"","")

library(lattice)
library(gplots)


pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Suppl_Figure4.pdf", height=8, width=5)
levelplot(z~x*y,grid, border="blue",cex.main=0.5, col.regions=rev(colorRampPalette(c("red","green", "aliceblue"), bias=0.2,interpolate="linear")(200)),xlab="",ylab="",at=c(0,1.30,seq(2,12,0.7)), colorkey=list(space="right"),scales=list(x=list(at=seq(1:ncol(y)),labels=colnames(y),rot=45),y = list(at=seq(1:nrow(y)),labels=rownames(y),font=3)),
panel=function(...) {
                       arg <- list(...)
                       panel.levelplot(...,)
                       panel.text(grid$x, grid$y,grid$sign)})
dev.off()






####################################################################
##### PERMUTATION TO SEE IF ENRICHED P-VALUE IN 46 SNPs FOR CHD ####
####################################################################

### READ A RESULTS OF META-ANALYSIS, JUST TO HAVE ALL THE SNPs in the ANALYSIS

ph <- "MG18_2"
cc <- read.table(paste("/proj/b2011036/nobackup/andrea/CHD_metabo/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/andrea/CHD_metabo/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table(paste("/proj/b2011036/nobackup/andrea/CHD_metabo/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)

ref <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/snp138_mod", sep=" ", stringsAsFactor=F)
ref$pos <- paste(ref[,1],":",(as.numeric(ref[,2])+1),sep="")

#### MERGE WITH THE SNPid ###
tm <- merge(t,ref,by.x="MarkerName",by.y="pos")

#### MERGE WITH HAPMAP SNPS ###
hap <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/hapmap_CEU_r23a.bim", sep="\t")

write.table(data.frame(SNP=tm$V3[tm$V3%in%hap$V2],P=runif(length(tm$V3[tm$V3%in%hap$V2]))), file="/proj/b2011036/nobackup/andrea/CHD_metabo/snp_for_pruning_all", quote=F, row.names=F)



#### OBTAINED PRUNED SNP DATASET TO EXTRACT 46 RANDOM SNPs ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--clump snp_for_pruning_all \
--clump-r2 0.2 --clump-p1 1 --clump-p2 1  \
--out snp_pruned_all


### READ PRUNED SNPs ###
prun <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/snp_pruned_all.clumped", header=T, stringsAsFactor=F)

### Read 44 CHD SNPs + P-values ###
p_noadj <- read.csv("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/fig4_noadj.csv", stringsAsFactor=F)
p_noadj <- apply(p_noadj[!p_noadj$X%in%c("PLA2G6","LPL","FADS1/FADS2","LIPG","LIPA","LCAT"),2:5],2,function(x){10^-as.numeric(x)})
props_noadj <- apply(p_noadj,2,function(x){sum(x<0.05)})

p_adj <- read.csv("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/fig4_adj.csv", stringsAsFactor=F)
p_adj <- apply(p_adj[!p_adj$X%in%c("PLA2G6","LPL","FADS1/FADS2","LIPG","LIPA","LCAT"),2:5],2,function(x){10^-as.numeric(x)})
props_adj <- apply(p_adj,2,function(x){sum(x<0.05)})


### MG 18:2 ###

### Keep only pruned SNPs ###
tmp <- tm[tm$V3%in%prun[,3],]
tmp <- tmp[!tmp$V3%in%ann$V2,]

### Hypergeometric test ###
phyper(props_adj[3],sum(tmp$P.value<0.05),nrow(tmp)-sum(tmp$P.value<0.05),44,lower.tail = F)
phyper(props_noadj[3],sum(tmp$P.value<0.05),nrow(tmp)-sum(tmp$P.value<0.05),44,lower.tail = F)


### LPC 18:2 ####

ph <- "LPC18_2"
cc <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)

#### MERGE WITH THE SNPid ###
tm2 <- merge(t,ref,by.x="MarkerName",by.y="pos")

### Keep only pruned SNPs ###
tmp2 <- tm2[tm2$V3%in%prun[,3],]
tmp2 <- tmp2[!tmp2$V3%in%ann$V2,]
### Hypergeometric test ###
phyper(props_adj[2],sum(tmp2$P.value<0.05),nrow(tmp2)-sum(tmp2$P.value<0.05),44,lower.tail = F)
phyper(props_noadj[2],sum(tmp2$P.value<0.05),nrow(tmp2)-sum(tmp2$P.value<0.05),44,lower.tail = F)



### LPC 18:1 ####

ph <- "LPC18_1"
cc <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)

#### MERGE WITH THE SNPid ###
tm3 <- merge(t,ref,by.x="MarkerName",by.y="pos")

### Keep only pruned SNPs ###
tmp3 <- tm3[tm3$V3%in%prun[,3],]
tmp3 <- tmp3[!tmp3$V3%in%ann$V2,]
### Hypergeometric test ###
phyper(props_adj[1],sum(tmp3$P.value<0.05),nrow(tmp3)-sum(tmp3$P.value<0.05),44,lower.tail = F)
phyper(props_adj[1],sum(tmp3$P.value<0.05),nrow(tmp3)-sum(tmp3$P.value<0.05),44,lower.tail = F)


### PE-cer ####

ph <- "PE_cer"
cc <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)

#### MERGE WITH THE SNPid ###
tm4 <- merge(t,ref,by.x="MarkerName",by.y="pos")

### Keep only pruned SNPs ###
tmp4 <- tm4[tm4$V3%in%prun[,3],]
tmp4 <- tmp4[!tmp4$V3%in%ann$V2,]
### Hypergeometric test ###
phyper(props_adj[4],sum(tmp4$P.value<0.05),nrow(tmp4)-sum(tmp4$P.value<0.05),44,lower.tail = F)
phyper(props_adj[4],sum(tmp4$P.value<0.05),nrow(tmp4)-sum(tmp4$P.value<0.05),44,lower.tail = F)





###########################################################################
##### NOT INCLUDED IN THE PAPER - RARE VARIANTS EXOME CHIP ANALYSIS #######
###########################################################################




### Run in Bash to extract mutations for the gene list ###
phs=(LIPA LCAT LIPC LIPG LPL PLA2G6 MGLL LIPE PLA2G1B PLA2G10 FADS1 FADS2)
rm /home/andrea/glob/metabo_chd/exome_temp/list_selected_exom
for ph in ${phs[@]} ; do
for chr in `seq 1 22` ;do
sed '1,3d' /home/andrea/glob/metabo_chd/exome_temp/ESP6500SI-V2-SSA137.updatedProteinHgvs.chr${chr}.snps_indels.txt | grep ${ph} | awk '{print $1, $2, $13, $15, $22}' >> /home/andrea/glob/metabo_chd/exome_temp/list_selected_exom
done
done


### in R ###


ph <- c("LIPA","LCAT","LIPC","LIPG","LPL","PLA2G6","MGLL","LIPE","PLA2G1B","PLA2G10","FADS1","FADS2")

## Read list of mutations 
d <- read.table("/home/andrea/glob/metabo_chd/exome_temp/list_selected_exom", stringsAsFactor=F)
d$pos <- sapply(strsplit(d[,1],":"),"[[",2)

## Keep only mutations with exactely the same label and missense only 
ds <- d[d[,3]%in%ph & d[,4]=="missense",]

## Match for position with the exome chip .bim file the position
a <- read.table("/home/andrea/glob/metabo_chd/exome_temp/exomeChip-ulsam-pivus.zcall.CLEAN.bim", stringsAsFactor=F)

## Write SNPs to extract
write.table(a[a[,4]%in%ds$pos,2],file="/home/andrea/glob/metabo_chd/exome_temp/to_extract.txt", row.names=F,quote=F,col.names=F)

## Run in plink the SNP extraction
plink --noweb --bfile /lynx/cvol/v38/b2011036/pivus_ulsam.exome_chip/called/exomeChip-ulsam-pivus.zcall.CLEAN --extract /home/andrea/glob/metabo_chd/exome_temp/to_extract.txt --recodeA --out /home/andrea/glob/metabo_chd/exome_temp/extracted_exome

### in R again ####

## Read the extracted SNPs
ext <- read.table("/home/andrea/glob/metabo_chd/exome_temp/extracted_exome.raw", header=T)


### Delete those with less than 10 subjects ###
exts <- ext[,7:ncol(ext)]
extss <- exts[,colSums(exts, na.rm=T)>10]
extsss <- cbind(ext[,1],extss)
colnames(extsss)[1] <- "FID"


## Merge the two studies
mULSAMPIVUS <- rbind(mULSAMF[,c("gwas_id","LPC18_1","LPC18_2","MG18_2","LPC18_2.hdl")],mPIVUSF[,c("gwas_id","LPC18_1","LPC18_2","MG18_2","LPC18_2.hdl")])
mULSAMPIVUS2 <- merge(extsss,mULSAMPIVUS,by.x="FID", by.y="gwas_id")


## Run association analysis
PRES181 <- NULL
PRES182 <- NULL
PRESmg <- NULL
PRES182hdl <- NULL

for(i in colnames(extsss)[2:ncol(extsss)])
{
	temp181 <- lm(scale(as.numeric(mULSAMPIVUS2$LPC18_1))~mULSAMPIVUS2[,i] ,data=mULSAMPIVUS2)
	temp182 <- lm(scale(as.numeric(mULSAMPIVUS2$LPC18_2))~mULSAMPIVUS2[,i] ,data=mULSAMPIVUS2)
	tempmg <-  lm(scale(as.numeric(mULSAMPIVUS2$MG18_2))~mULSAMPIVUS2[,i] ,data=mULSAMPIVUS2)
	temp182hdl <- lm(scale(as.numeric(mULSAMPIVUS2$LPC18_2.hdl))~mULSAMPIVUS2[,i],data=mULSAMPIVUS2)
	
	PRES181 <- rbind(PRES181,c(i,temp181$coefficients[2],summary(temp181)$coefficients[2,2],summary(temp181)$coefficients[2,4],as.numeric(table(mULSAMPIVUS2[,i])[1]),as.numeric(table(mULSAMPIVUS2[,i])[2]),as.numeric(table(mULSAMPIVUS2[,i])[3])))
	PRES182 <- rbind(PRES182,c(i,temp182$coefficients[2],summary(temp182)$coefficients[2,2],summary(temp182)$coefficients[2,4]))
	PRESmg <- rbind(PRESmg,c(i,tempmg$coefficients[2],summary(tempmg)$coefficients[2,2],summary(tempmg)$coefficients[2,4]))
	PRES182hdl <- rbind(PRES182hdl,c(i,temp182hdl$coefficients[2],summary(temp182hdl)$coefficients[2,2],summary(temp182hdl)$coefficients[2,4]))
}


colnames(PRES181) <- c("snp","b181","se181","p181","0","1","2")
colnames(PRES182) <- c("snp","b182","se182","p182")
colnames(PRESmg) <- c("snp","bmg","semg","pmg")
colnames(PRES182hdl) <- c("snp","b182hdl","se182hdl","p182hdl")


res <- data.frame(cbind(sapply(strsplit(PRES181[,1],"_"),"[[",1), PRES181[,1:4],PRES182[,2:4],PRESmg[,2:4],PRES182hdl[,2:4],PRES181[,5:7]))







### Match results with the annotation file
a <- read.table("/home/andrea/glob/metabo_chd/exome_temp/exomeChip-ulsam-pivus.zcall.CLEAN.bim", stringsAsFactor=F)
d <- read.table("/home/andrea/glob/metabo_chd/exome_temp/list_selected_exom", stringsAsFactor=F)
d$pos <- sapply(strsplit(d[,1],":"),"[[",2)

resfin <- merge(res,a,by.x="V1",by.y="V2")
resfin2 <- merge(resfin,d,by.x="V4",by.y="pos")

write.csv(resfin2,"test.csv")



