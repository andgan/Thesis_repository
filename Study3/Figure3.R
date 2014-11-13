#----------------------------------------------
# Filename: Figure3.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 31MAR2014
# Updated: 
# Purpose: Mendelian randomization analysis using gtx package
# Note: 
#-----------------------------------------------
# Data used: LPC18_1_metaanalysis1.tbl LPC18_2_metaanalysis1.tbl MG18_2_metaanalysis1.tbl PE_cer_metaanalysis1.tbl snp138_mod CARDIoGRAM_GWAS_RESULTS.txt gwas_catalog.txt
# Data created: LPC_18_1_00005.pdf LPC_18_2_00005.pdf MG18_2_00005.pdf PECER_00005.pdf
#    LPC_18_1_check.pdf LPC_18_2_check.pdf MG18_2_check.pdf PECER_check.pdf
#-----------------------------------------------
# OP: R 2.13.1, plink
#-----------------------------------------------*/



library(gtx)		

#######################
### Lyso PC: 18_1 #####
#######################


	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/LPC18_1_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table("/proj/b2011036/nobackup/LPC18_1_metaanalysis1.tbl", header=T, comment.char = "", nrow=as.numeric(nrows), stringsAsFactor=F)


### READ TABLE 1000G SNPs b37 ####
ref <- read.table("/proj/b2011036/nobackup/snp138_mod", sep=" ", stringsAsFactor=F)
ref$pos <- paste(ref[,1],":",(as.numeric(ref[,2])+1),sep="")

#### MERGE WITH THE SNPid ###
tm <- merge(t,ref,by.x="MarkerName",by.y="pos")

### READ CARDIOGRAM RESULTS ####
cardio <- read.table("/home/andrea/glob/metabo_chd/CARDIoGRAM_GWAS_RESULTS.txt", header=T, stringsAsFactor=F)

#### Only SNPs with P < 5x10-5 AND ONLY SNPS in CARDIOGRAM ###
tms2 <- tm[as.numeric(tm$P.value) < 0.00005 & tm$V3%in%cardio$SNP,]


#### SAVE SNP FOR PRUNING ####
write.table(data.frame(SNP=tms2$V3,P=tms2$P.value), file="/proj/b2011036/nobackup/snp_for_pruning_LPC1812", quote=F, row.names=F)


#### RUN IN PLINK ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--clump snp_for_pruning_LPC1812 \
--clump-r2 0.2 --clump-p1 1 --clump-p2 1  \
--out snp_pruned_LPC1812
	

### READ PRUNED SNPs ###
prun2 <- read.table("/proj/b2011036/nobackup/snp_pruned_LPC1812.clumped", header=T, stringsAsFactor=F)


#### Cardio which is also in clumped ###
cardioP2 <- cardio[cardio$SNP%in%prun2$SNP, ]
tmP2 <- tm[tm$V3 %in% prun2$SNP, ]

### Now sort so that is the same of cardio ###
tmPS2 <- tmP2[order(match(tmP2$V3,cardioP2$SNP)),]

### Check risk alelle are the same ###
sum(toupper(as.character(tmPS2$Allele1))!=cardioP2$reference_allele & toupper(as.character(tmPS2$Allele1))!=cardioP2$other_allele)

### Now code that the risk allele is the same ###
tmPS2$beta_cor <- ifelse(toupper(as.character(tmPS2$Allele1))==cardioP2$reference_allele,tmPS2$Effect,-tmPS2$Effect)

### Save SNPs so that I can find proxy using SNAP ###
write.table(cardioP2$SNP,"/proj/b2011036/nobackup/for_snap_lpc181", row.names=F,quote=F,col.names=F)


#### RUN IN PLINK ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--r2 \
--ld-snp-list for_snap_lpc181 \
--ld-window-kb 5000 \
--ld-window-r2 0.2 \
--out proxy_LPC181



### READ GWAS catalog ###
gwas_c <- read.delim("/proj/b2011036/nobackup/gwas_catalog.txt", header=T, stringsAsFactor=F)

### READ results from Plink ###
proxies <- read.table("/proj/b2011036/nobackup/proxy_LPC181.ld", header=T, stringsAsFactor=F)


## Check traits to exclude
traits_to_exclude <- gwas_c$Disease.Trait[(gwas_c$SNPs%in%cardioP2$SNP) | (gwas_c$SNPs%in%proxies$SNP_B)]


## Check gwas_catalog for control
write.csv(gwas_c[gwas_c$SNPs%in%proxies$SNP_B,],file="test.csv")

traits_to_exclude <- traits_to_exclude[!traits_to_exclude%in%c("Migraine without aura","Age-related macular degeneration","Glaucoma")]

### Select SNPs to exclude because in GWAS catalog
gwas_c2 <- gwas_c[gwas_c$Disease.Trait%in%traits_to_exclude,]
to_excl_proxy <- gwas_c2$SNPs[gwas_c2$SNPs%in%proxies$SNP_B]
to_excl <- unique(proxies$SNP_A[proxies$SNP_B%in%to_excl_proxy])


## Exclusion based on heterogenity
yQ1 <- grs.filter.Qrs(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se)
## Exclusion based on GWAS catalog
## Additionally esclud rs7035484 because outlier
yQ2 <- ifelse(cardioP2$SNP%in%c(to_excl,"rs7035484"),F,T)


### Check if there is any outlier ###
IVestimate <- cardioP2$log_odds/tmPS2$beta_cor
seIVestimate <- abs(IVestimate)*sqrt((tmPS2$StdErr/tmPS2$beta_cor)^2+(cardioP2$log_odds_se/cardioP2$log_odds)^2)

col <- ifelse(cardioP2$SNP%in%c(to_excl,"rs7035484"),"red","black")
library(gplots)
pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/LPC_18_1_check.pdf")
plotCI(IVestimate,ui=IVestimate+1.96*seIVestimate,li=IVestimate-1.96*seIVestimate, ylab="IV estimates",xaxt = "n", xlab="", col=col, main="LysoPC 18:1")
abline(h=0)
axis(1,at=1:length(IVestimate),labels=cardioP2$SNP,las=2, cex.axis=0.8)
dev.off()


y <- grs.summary(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se, cardioP2$N_case+cardioP2$N_control)
y1 <- grs.summary(tmPS2$beta_cor[yQ1], cardioP2$log_odds[yQ1], cardioP2$log_odds_se[yQ1], cardioP2$N_case[yQ1]+cardioP2$N_control[yQ1])
y2 <- grs.summary(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$N_case[yQ2]+cardioP2$N_control[yQ2])



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/LPC_18_1_00005.pdf", height=7,width=7)
grs.plot(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$SNP[yQ2],textcex=0.8)
title(ylab = "ln(odds) change in CHD risk score per allele", xlab = "ln change in 1-SD of LysoPC 18:1 per allele", main="LysoPC 18:1", cex.main=1.4,cex.lab=1)
dev.off()





#####################
### LysoPC 18:2 #####
#####################
	

	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/LPC18_2_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table("/proj/b2011036/nobackup/LPC18_2_metaanalysis1.tbl", header=T, comment.char = "", nrow=as.numeric(nrows), stringsAsFactor=F)


### READ TABLE 1000G SNPs b37 ####
ref <- read.table("/proj/b2011036/nobackup/snp138_mod", sep=" ", stringsAsFactor=F)
ref$pos <- paste(ref[,1],":",(as.numeric(ref[,2])+1),sep="")

#### MERGE WITH THE SNPid ###
tm <- merge(t,ref,by.x="MarkerName",by.y="pos")

### READ CARDIOGRAM RESULTS ####
cardio <- read.table("/home/andrea/glob/metabo_chd/CARDIoGRAM_GWAS_RESULTS.txt", header=T, stringsAsFactor=F)

#### Only SNPs with P < 5x10-5 AND ONLY SNPS in CARDIOGRAM ###
tms2 <- tm[as.numeric(tm$P.value) < 0.00005 & tm$V3%in%cardio$SNP,]


#### SAVE SNP FOR PRUNING ####
write.table(data.frame(SNP=tms2$V3,P=tms2$P.value), file="/proj/b2011036/nobackup/snp_for_pruning_LPC1822", quote=F, row.names=F)


#### RUN IN PLINK ####

plink --bfile hapmap_CEU_r23a \
--noweb \
--clump snp_for_pruning_LPC1822 \
--clump-r2 0.2 --clump-p1 1 --clump-p2 1  \
--out snp_pruned_LPC1822


### READ PRUNED SNPs ###
prun2 <- read.table("/proj/b2011036/nobackup/snp_pruned_LPC1822.clumped", header=T, stringsAsFactor=F)

#### Cardio which is also in clumped ###
cardioP2 <- cardio[cardio$SNP%in%prun2$SNP, ]
tmP2 <- tm[tm$V3 %in% prun2$SNP, ]

### Now sort so that is the same of cardio ###
tmPS2 <- tmP2[order(match(tmP2$V3,cardioP2$SNP)),]

### Check risk alelle are the same ###
sum(toupper(as.character(tmPS2$Allele1))!=cardioP2$reference_allele & toupper(as.character(tmPS2$Allele1))!=cardioP2$other_allele)

### Now code that the risk allele is the same ###
tmPS2$beta_cor <- ifelse(toupper(as.character(tmPS2$Allele1))==cardioP2$reference_allele,tmPS2$Effect,-tmPS2$Effect)

### Save SNPs so that I can find proxy using SNAP ###
write.table(cardioP2$SNP,"/proj/b2011036/nobackup/for_snap_lpc182", row.names=F,quote=F,col.names=F)


#### RUN IN PLINK ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--r2 \
--ld-snp-list for_snap_lpc182 \
--ld-window-kb 5000 \
--ld-window-r2 0.2 \
--out proxy_LPC182



### READ GWAS catalog ###
gwas_c <- read.delim("/proj/b2011036/nobackup/gwas_catalog.txt", header=T, stringsAsFactor=F)

### READ results from Plink ###
proxies <- read.table("/proj/b2011036/nobackup/proxy_LPC182.ld", header=T, stringsAsFactor=F)


## Check traits to exclude
traits_to_exclude <- gwas_c$Disease.Trait[(gwas_c$SNPs%in%cardioP2$SNP) | (gwas_c$SNPs%in%proxies$SNP_B)]

## Check gwas_catalog for control
write.csv(gwas_c[gwas_c$SNPs%in%proxies$SNP_B,],file="test.csv")

## This traits should be included because non cardio-metabolic traits
traits_to_exclude <- traits_to_exclude[!traits_to_exclude%in%c("Phospholipid levels (plasma)","Celiac disease","Age-related macular degeneration","Major depressive disorder","Glaucoma")]


### Select SNPs to exclude because in GWAS catalog
gwas_c2 <- gwas_c[gwas_c$Disease.Trait%in%traits_to_exclude,]
to_excl_proxy <- gwas_c2$SNPs[gwas_c2$SNPs%in%proxies$SNP_B]
to_excl <- unique(proxies$SNP_A[proxies$SNP_B%in%to_excl_proxy])


## Exclusion based on heterogenity
yQ1 <- grs.filter.Qrs(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se)
## Exclusion based on GWAS catalog
## We additionally exclude rs2106120 because outlier
yQ2 <- ifelse(cardioP2$SNP%in%c(to_excl,"rs2106120"),F,T)


### Check if there is any outlier ###
IVestimate <- cardioP2$log_odds/tmPS2$beta_cor
seIVestimate <- abs(IVestimate)*sqrt((tmPS2$StdErr/tmPS2$beta_cor)^2+(cardioP2$log_odds_se/cardioP2$log_odds)^2)

col <- ifelse(cardioP2$SNP%in%c(to_excl,"rs2106120"),"red","black")
library(gplots)
pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/LPC_18_2_check.pdf")
plotCI(IVestimate,ui=IVestimate+1.96*seIVestimate,li=IVestimate-1.96*seIVestimate, ylab="IV estimates",xaxt = "n", xlab="", col=col, main="LysoPC 18:2")
abline(h=0)
axis(1,at=1:length(IVestimate),labels=cardioP2$SNP,las=2, cex.axis=0.8)
dev.off()


y <- grs.summary(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se, cardioP2$N_case+cardioP2$N_control)
y1 <- grs.summary(tmPS2$beta_cor[yQ1], cardioP2$log_odds[yQ1], cardioP2$log_odds_se[yQ1], cardioP2$N_case[yQ1]+cardioP2$N_control[yQ1])
y2 <- grs.summary(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$N_case[yQ2]+cardioP2$N_control[yQ2])



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/LPC_18_2_00005.pdf", height=7,width=7)
grs.plot(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$SNP[yQ2],textcex=0.8)
title(ylab = "ln(odds) change in CHD risk score per allele", xlab = "ln change in 1-SD of LysoPC 18:2 per allele", main="LysoPC 18:2", cex.main=1.4,cex.lab=1)
dev.off()



#################
### MG 18_2 #####
#################


library(gtx)		

	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/MG18_2_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table("/proj/b2011036/nobackup/MG18_2_metaanalysis1.tbl", header=T, comment.char = "", nrow=as.numeric(nrows),  stringsAsFactor=F)


### READ TABLE 1000G SNPs b37 ####
ref <- read.table("/proj/b2011036/nobackup/snp138_mod", sep=" ", stringsAsFactor=F)
ref$pos <- paste(ref[,1],":",(as.numeric(ref[,2])+1),sep="")

#### MERGE WITH THE SNPid ###
tm <- merge(t,ref,by.x="MarkerName",by.y="pos")

### READ CARDIOGRAM RESULTS ####
cardio <- read.table("/home/andrea/glob/metabo_chd/CARDIoGRAM_GWAS_RESULTS.txt", header=T, stringsAsFactor=F)

#### Only SNPs with P < 5x10-5 AND ONLY SNPS in CARDIOGRAM ###
tms2<- tm[as.numeric(tm$P.value) < 0.00005 & tm$V3%in%cardio$SNP,]

#p.adjust(as.numeric(tm$P.value),method="BH")[tm$P.value==5.005e-05]
#tm$P.value[which(abs(tm$P.value-0.00005)==min(abs(tm$P.value-0.00005)))]


#### SAVE SNP FOR PRUNING ####
write.table(data.frame(SNP=tms2$V3,P=tms2$P.value), file="/proj/b2011036/nobackup/snp_for_pruning_MG18_22", quote=F, row.names=F)

#### RUN IN PLINK ####


plink --bfile hapmap_CEU_r23a \
--noweb \
--clump snp_for_pruning_MG18_22 \
--clump-r2 0.2 --clump-p1 1 --clump-p2 1  \
--out snp_pruned_MG18_22


	
### READ PRUNED SNPs ###
prun2 <- read.table("/proj/b2011036/nobackup/snp_pruned_MG18_22.clumped", header=T, stringsAsFactor=F)


#### Cardio which is also in clumped ###
cardioP2 <- cardio[cardio$SNP%in%prun2$SNP, ]
tmP2 <- tm[tm$V3 %in% prun2$SNP, ]

### Now sort so that is the same of cardio ###
tmPS2 <- tmP2[order(match(tmP2$V3,cardioP2$SNP)),]

### Check risk alelle are the same ###
sum(toupper(as.character(tmPS2$Allele1))!=cardioP2$reference_allele & toupper(as.character(tmPS2$Allele1))!=cardioP2$other_allele)

### Now code that the risk allele is the same ###
tmPS2$beta_cor <- ifelse(toupper(as.character(tmPS2$Allele1))==cardioP2$reference_allele,tmPS2$Effect,-tmPS2$Effect)

### Save SNPs so that I can find proxy using SNAP ###
write.table(cardioP2$SNP,"/proj/b2011036/nobackup/for_snap_mg182", row.names=F,quote=F,col.names=F)


#### RUN IN PLINK ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--r2 \
--ld-snp-list for_snap_mg182 \
--ld-window-kb 5000 \
--ld-window-r2 0.2 \
--out proxy_MG182



### READ GWAS catalog ###
gwas_c <- read.delim("/proj/b2011036/nobackup/gwas_catalog.txt", header=T, stringsAsFactor=F)

### READ results from Plink ###
proxies <- read.table("/proj/b2011036/nobackup/proxy_MG182.ld", header=T, stringsAsFactor=F)


## Check traits to exclude
traits_to_exclude <- gwas_c$Disease.Trait[(gwas_c$SNPs%in%cardioP2$SNP) | (gwas_c$SNPs%in%proxies$SNP_B)]

## Check gwas_catalog for control
write.csv(gwas_c[gwas_c$SNPs%in%proxies$SNP_B,],file="test.csv")

## This traits should be included because non cardio-metabolic traits
traits_to_exclude <- traits_to_exclude[!traits_to_exclude%in%c("Metabolite levels","Response to Vitamin E supplementation","Vitamin E levels","Coronary heart disease")]



### Select SNPs to exclude because in GWAS catalog
gwas_c2 <- gwas_c[gwas_c$Disease.Trait%in%traits_to_exclude,]
to_excl_proxy <- gwas_c2$SNPs[gwas_c2$SNPs%in%proxies$SNP_B]
to_excl <- unique(proxies$SNP_A[proxies$SNP_B%in%to_excl_proxy])


## Exclusion based on heterogenity
yQ1 <- grs.filter.Qrs(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se)
## Exclusion based on GWAS catalog
## We additionally exclude rs2106120 because outlier
yQ2 <- ifelse(cardioP2$SNP%in%to_excl,F,T)


### Check if there is any outlier ###
IVestimate <- cardioP2$log_odds/tmPS2$beta_cor
seIVestimate <- abs(IVestimate)*sqrt((tmPS2$StdErr/tmPS2$beta_cor)^2+(cardioP2$log_odds_se/cardioP2$log_odds)^2)

col <- ifelse(cardioP2$SNP%in%to_excl,"red","black")
library(gplots)
pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/MG_18_2_check.pdf")
plotCI(IVestimate,ui=IVestimate+1.96*seIVestimate,li=IVestimate-1.96*seIVestimate, ylab="IV estimates",xaxt = "n", xlab="", col=col, main="MG 18:2")
abline(h=0)
axis(1,at=1:length(IVestimate),labels=cardioP2$SNP,las=2, cex.axis=0.8)
dev.off()


y <- grs.summary(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se, cardioP2$N_case+cardioP2$N_control)
y1 <- grs.summary(tmPS2$beta_cor[yQ1], cardioP2$log_odds[yQ1], cardioP2$log_odds_se[yQ1], cardioP2$N_case[yQ1]+cardioP2$N_control[yQ1])
y2 <- grs.summary(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$N_case[yQ2]+cardioP2$N_control[yQ2])



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/MG_18_2_00005.pdf", height=7,width=7)
grs.plot(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$SNP[yQ2],textcex=0.8)
title(ylab = "ln(odds) change in CHD risk score per allele", xlab = "ln change in 1-SD of MG 18:2 per allele", main="MG 18:2", cex.main=1.4,cex.lab=1)
dev.off()




#################
### SM 28:1 #####
#################

library(gtx)		

	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/andrea/CHD_metabo/PE_cer_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/PE_cer_metaanalysis1.tbl", header=T, comment.char = "", nrow=as.numeric(nrows),  stringsAsFactor=F)


### READ TABLE 1000G SNPs b37 ####
ref <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/snp138_mod", sep=" ", stringsAsFactor=F)
ref$pos <- paste(ref[,1],":",(as.numeric(ref[,2])+1),sep="")

#### MERGE WITH THE SNPid ###
tm <- merge(t,ref,by.x="MarkerName",by.y="pos")

### READ CARDIOGRAM RESULTS ####
cardio <- read.table("/home/andrea/glob/metabo_chd/CARDIoGRAM_GWAS_RESULTS.txt", header=T, stringsAsFactor=F)

#### Only SNPs with P < 5x10-5 AND ONLY SNPS in CARDIOGRAM ###
tms2 <- tm[as.numeric(tm$P.value) < 0.00005 & tm$V3%in%cardio$SNP,]

#### SAVE SNP FOR PRUNING ####
write.table(data.frame(SNP=tms2$V3,P=tms2$P.value), file="/proj/b2011036/nobackup/andrea/CHD_metabo/snp_for_pruning_PE_cer2", quote=F, row.names=F)


#### RUN IN PLINK ####
	
plink --bfile hapmap_CEU_r23a \
--noweb \
--clump snp_for_pruning_PE_cer2 \
--clump-r2 0.2 --clump-p1 1 --clump-p2 1  \
--out snp_pruned_PE_cer2

	

### READ PRUNED SNPs ###
prun2 <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/snp_pruned_PE_cer2.clumped", header=T, stringsAsFactor=F)


#### Cardio which is also in clumped ###
cardioP2 <- cardio[cardio$SNP%in%prun2$SNP, ]
tmP2 <- tm[tm$V3 %in% prun2$SNP, ]

### Now sort so that is the same of cardio ###
tmPS2 <- tmP2[order(match(tmP2$V3,cardioP2$SNP)),]

### Check risk alelle are the same ###
sum(toupper(as.character(tmPS2$Allele1))!=cardioP2$reference_allele & toupper(as.character(tmPS2$Allele1))!=cardioP2$other_allele)

### Now code that the risk allele is the same ###
tmPS2$beta_cor <- ifelse(toupper(as.character(tmPS2$Allele1))==cardioP2$reference_allele,tmPS2$Effect,-tmPS2$Effect)

### Save SNPs so that I can find proxy using SNAP ###
write.table(cardioP2$SNP,"/proj/b2011036/nobackup/andrea/CHD_metabo/for_snap_pe_cer", row.names=F,quote=F,col.names=F)


#### RUN IN PLINK ####
plink --bfile hapmap_CEU_r23a \
--noweb \
--r2 \
--ld-snp-list for_snap_pe_cer \
--ld-window-kb 5000 \
--ld-window-r2 0.2 \
--out proxy_PECER



### READ GWAS catalog ###
gwas_c <- read.delim("/proj/b2011036/nobackup/andrea/CHD_metabo/gwas_catalog.txt", header=T, stringsAsFactor=F)

### READ results from Plink ###
proxies <- read.table("/proj/b2011036/nobackup/andrea/CHD_metabo/proxy_PECER.ld", header=T, stringsAsFactor=F)


## Check traits to exclude
traits_to_exclude <- gwas_c$Disease.Trait[(gwas_c$SNPs%in%cardioP2$SNP) | (gwas_c$SNPs%in%proxies$SNP_B)]

## Check gwas_catalog for control
write.csv(gwas_c[gwas_c$SNPs%in%proxies$SNP_B,],file="test.csv")


### Select SNPs to exclude because in GWAS catalog
gwas_c2 <- gwas_c[gwas_c$Disease.Trait%in%traits_to_exclude,]
to_excl_proxy <- gwas_c2$SNPs[gwas_c2$SNPs%in%proxies$SNP_B]
to_excl <- unique(proxies$SNP_A[proxies$SNP_B%in%to_excl_proxy])


## Exclusion based on heterogenity
yQ1 <- grs.filter.Qrs(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se)
## Exclusion based on GWAS catalog
## We additionally exclude rs2106120 because outlier
yQ2 <- ifelse(cardioP2$SNP%in%to_excl,F,T)


### Check if there is any outlier ###
IVestimate <- cardioP2$log_odds/tmPS2$beta_cor
seIVestimate <- abs(IVestimate)*sqrt((tmPS2$StdErr/tmPS2$beta_cor)^2+(cardioP2$log_odds_se/cardioP2$log_odds)^2)

col <- ifelse(cardioP2$SNP%in%to_excl,"red","black")
library(gplots)
pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/PECER_check.pdf")
plotCI(IVestimate,ui=IVestimate+1.96*seIVestimate,li=IVestimate-1.96*seIVestimate, ylab="IV estimates",xaxt = "n", xlab="", col=col, main="SM 34:1")
abline(h=0)
axis(1,at=1:length(IVestimate),labels=cardioP2$SNP,las=2, cex.axis=0.8)
dev.off()


y <- grs.summary(tmPS2$beta_cor, cardioP2$log_odds, cardioP2$log_odds_se, cardioP2$N_case+cardioP2$N_control)
y1 <- grs.summary(tmPS2$beta_cor[yQ1], cardioP2$log_odds[yQ1], cardioP2$log_odds_se[yQ1], cardioP2$N_case[yQ1]+cardioP2$N_control[yQ1])
y2 <- grs.summary(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$N_case[yQ2]+cardioP2$N_control[yQ2])



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/PECER_00005.pdf", height=7,width=7)
grs.plot(tmPS2$beta_cor[yQ2], cardioP2$log_odds[yQ2], cardioP2$log_odds_se[yQ2], cardioP2$SNP[yQ2],textcex=0.8)
title(ylab = "ln(odds) change in CHD risk score per allele", xlab = "ln change in 1-SD of SM 28:1 per allele", main="SM 28:1", cex.main=1.4,cex.lab=1)
dev.off()

