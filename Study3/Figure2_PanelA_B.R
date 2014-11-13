#----------------------------------------------
# Filename: Figure2_PanelA_B.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 12JUL2014 - Updated name of the tables and covariates used for adjustment/
# Purpose: Analysis used in the paper. Figure 2, panel A and B
# Note: 
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small, Pivus small), PIVUS pek tander vs athero6.txt PIVUS data PLA2s.txt
# Data created: Fig_2_p1a.pdf Fig_2_p1b.pdf Fig_2_p1c.pdf Fig_2_p1d.pdf Fig_2_p2a.pdf Fig_2_p2b.pdf
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/


###############
#### INDEX ####
### 1. FIGURE 2, PANEL A - ASSOCIATION BETWEEN 4 metabolites and established risk factors 
### 2. FIGURE 2, PANEL B - ASSOCIATION BETWEEN 4 metabolites and inflammatory markers and subclinical CVD in PIVUS




#############################################################################################
## 1. FIGURE 2, PANEL A - ASSOCIATION BETWEEN 4 metabolites and established risk factors ####
#############################################################################################

## PIVUS ##
load("/home/andrea/glob/alignment_pivus_small/Results/Final_datasets/step4.Rdata")

## 

# M496.340T365.792 = LysoPC 16:0
# M524.372T425.563 = LysoPC 18:0
# M522.356T384.872 = LysoPC 18:1
# M520.340T347.575 = LysoPC 18:2
# M566.322T352.005 = LysoPC 20:4
# M377.267T388.196 = MG 18:2
# M661.528T579.938 = PE_cer



mPIVUS <- metabo_p[,c("id","M496.340T365.792","M524.372T425.563","M522.356T384.872","M520.340T347.575","M566.322T352.005",
"M377.267T388.196","M661.528T579.938","check_d","check_d","smoke01","diab","tg","antihyp")]
cardiphen <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS pek tander vs athero6.txt", header=T, sep="\t", stringsAsFactor=F)
pla2 <- read.table("/home/andrea/glob/alignment_pivus_small/Data/PIVUS data PLA2s.txt", header=T, sep="\t", stringsAsFactor=F)


# Numeric
mPIVUS$id <- as.numeric(mPIVUS$id)
pla2 <- apply(pla2,2,as.numeric)

# Merging
tmp <- merge(mPIVUS,cardiphen,by="id")
mPIVUSF <- merge(tmp,pla2,by.x="id", by.y="lpnr")

colnames(mPIVUSF)[2:8] <- c("LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC20_4","MG18_2","PE_cer")

mPIVUSF$logcrp <- log(mPIVUSF$crp)
mPIVUSF$logtg <- log(mPIVUSF$tg)

mPIVUSF$lnhomocyst <- log(mPIVUSF$homocyst)
mPIVUSF$Lngssg <- log(mPIVUSF$gssg)
mPIVUSF$Lngsh <- log(mPIVUSF$gsh)
mPIVUSF$gluratio <- mPIVUSF$Lngssg/mPIVUSF$Lngsh




## TWINGENE ##


load("/home/andrea/glob/alignment_twge_small/Results/Final_datasets/step4.Rdata")

## 
# M496.340T364.064 = LysoPC 16:0
# M524.372T423.532 = LysoPC 18:0
# M522.356T382.722 = LysoPC 18:1
# M520.340T346.671 = LysoPC 18:2
# M566.322T340.659 = LysoPC 20:4
# M393.231T384.578 = MG 18;2
# M661.528T580.694 = PE_cer

mTWINGENEF <- metabo_p[,c("twinnr","M496.340T364.064","M524.372T423.532","M522.356T382.722","M520.340T346.671","M566.322T340.659","M393.231T384.578","M661.528T580.694",colnames(metabo_p)[9756:ncol(metabo_p)])]

colnames(mTWINGENEF)[2:8] <- c("LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC20_4","MG18_2","PE_cer")

mTWINGENEF$logcrp <- log(mTWINGENEF$crp)
mTWINGENEF$logtg <- log(mTWINGENEF$tg)



## ULSAM ##
load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

## 
# M496.340T375.465 = LysoPC 16:0
# M524.372T436.589 = LysoPC 18:0
# M522.356T395.236 = LysoPC 18:1
# M520.340T357.723 = LysoPC 18:2
# M566.322T351.269 = LysoPC 20:4
# M393.238T397.936 = MG 18:2
# M661.528T590.990 = PE_cer



mULSAMF <- metabo_p[metabo_p$time==0 & metabo_p$double_==0,c("pat_time","M496.340T375.465","M524.372T436.589","M522.356T395.236","M520.340T357.723","M566.322T351.269","M393.238T397.936","M661.528T590.990",colnames(metabo_p)[10163:ncol(metabo_p)])]

colnames(mULSAMF)[2:8] <- c("LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC20_4","MG18_2","PE_cer")

mULSAMF$logcrp <- log(mULSAMF$crp)
mULSAMF$logtg <- log(mULSAMF$tg)




### FHS RISK FACTORS (ULSAM) ###

lysoPCL <- c("LPC18_1","LPC18_2","MG18_2","PE_cer")
## Association
resfhsU <- list(NA)
for (k in lysoPCL)
{
		temp <- NULL
		for (i in rev(c("ldl","hdl","logtg","bmi","sbp")))
			{
				st_r <- scale(as.numeric(mULSAMF[,i]))
				mod <- lm(st_r~scale(as.numeric(mULSAMF[,k])))
				temp <- rbind(temp,c(i,coefficients(mod)[2],coefficients(summary(mod))[2,2],coefficients(summary(mod))[2,4]))
			}
	resfhsU[[which(k==lysoPCL)]] <- temp
}




### FHS RISK FACTORS (TWGE) ###

lysoPCL <- c("LPC18_1","LPC18_2","MG18_2","PE_cer")
## Association
resfhsT <- list(NA)
for (k in lysoPCL)
{
temp <- NULL
	for (i in rev(c("ldl","hdl","logtg","bmi","sbp")))
		{
			st_r <- scale(as.numeric(mTWINGENEF[,i]))
			mod <- lm(st_r~mTWINGENEF$age+scale(as.numeric(mTWINGENEF[,k]))+mTWINGENEF$sex)
			temp <- rbind(temp,c(i,coefficients(mod)[3],coefficients(summary(mod))[3,2],coefficients(summary(mod))[3,4]))
		}
	resfhsT[[which(k==lysoPCL)]] <- temp
}


### FHS RISK FACTORS (PIVUS) ###

lysoPCL <- c("LPC18_1","LPC18_2","MG18_2","PE_cer")
## Association
resfhsP <- list(NA)
for (k in lysoPCL)
{
temp <- NULL
	for (i in rev(c("ldl","hdl","logtg","bmi","manuelltsbp")))
		{
			st_r <- scale(as.numeric(mPIVUSF[,i]))
			mod <- lm(st_r~mPIVUSF$kvinna1+scale(as.numeric(mPIVUSF[,k])))
			temp <- rbind(temp,c(i,coefficients(mod)[3],coefficients(summary(mod))[3,2],coefficients(summary(mod))[3,4]))
		}
	resfhsP[[which(k==lysoPCL)]] <- temp
}


### Established risk factors ###

fhsuadj <- NULL
for (i in 1:nrow(resfhsU[[1]]))
{
	fhsuadj <- rbind(fhsuadj,rbind(resfhsU[[1]][i,],resfhsT[[1]][i,],resfhsP[[1]][i,]))
}

fhsuadj2 <- data.frame(as.numeric(fhsuadj[,2]),as.numeric(fhsuadj[,2])-1.96*as.numeric(fhsuadj[,3]),as.numeric(fhsuadj[,2])+1.96*as.numeric(fhsuadj[,3]))
fhsuadj2$x <- c(1,2,3,6,7,8,11,12,13,16,17,18,21,22,23)
colnames(fhsuadj2) <- c("y","ylo","yhi","x")



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p1a.pdf", height=6, width=5)
ggplot(fhsuadj2, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=c(rep(c("red","blue","green"),5))))+
geom_pointrange(size=1, shape=22)+
coord_flip() + 
geom_hline(aes(x=0), lty=2) + 
labs(x = "", y = "Beta for 1 SD increase in LysoPC 18:1") + theme(axis.text.y=element_text(colour="black",size = rel(1.2)),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size = rel(1)),legend.position=c(.8, .15)) + 
ylim(-0.6,+0.6) + 
ggtitle("LysoPC 18:1") + 
scale_x_discrete(breaks=c(2,7,12,17,22),labels=rev(c("LDL-cholesterol","HDL-cholesterol","Log Triglycerides","BMI","Systolic blood pressure")))+
scale_fill_discrete(breaks=rev(c("red","blue","green")),labels=c("Ulsam","TwinGene","PIVUS"))+
guides(fill=guide_legend(title=NULL)) 
dev.off()


fhsuadj <- NULL
for (i in 1:nrow(resfhsU[[2]]))
{
	fhsuadj <- rbind(fhsuadj,rbind(resfhsU[[2]][i,],resfhsT[[2]][i,],resfhsP[[2]][i,]))
}

fhsuadj2 <- data.frame(as.numeric(fhsuadj[,2]),as.numeric(fhsuadj[,2])-1.96*as.numeric(fhsuadj[,3]),as.numeric(fhsuadj[,2])+1.96*as.numeric(fhsuadj[,3]))
fhsuadj2$x <- c(1,2,3,6,7,8,11,12,13,16,17,18,21,22,23)
colnames(fhsuadj2) <- c("y","ylo","yhi","x")



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p1b.pdf", height=6, width=4)
ggplot(fhsuadj2, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=c(rep(c("red","blue","green"),5))))+
geom_pointrange(size=1, shape=22)+
coord_flip() + 
geom_hline(aes(x=0), lty=2) + 
labs(x = "", y = "Beta for 1 SD increase in LysoPC 18:2") + theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size = rel(1))) + guides(fill=FALSE) +
ylim(-0.6,+0.6) + 
ggtitle("LysoPC 18:2") +
scale_x_discrete(breaks=c(2,7,12,17,22))
dev.off()




fhsuadj <- NULL
for (i in 1:nrow(resfhsU[[3]]))
{
	fhsuadj <- rbind(fhsuadj,rbind(resfhsU[[3]][i,],resfhsT[[3]][i,],resfhsP[[3]][i,]))
}

fhsuadj2 <- data.frame(as.numeric(fhsuadj[,2]),as.numeric(fhsuadj[,2])-1.96*as.numeric(fhsuadj[,3]),as.numeric(fhsuadj[,2])+1.96*as.numeric(fhsuadj[,3]))
fhsuadj2$x <- c(1,2,3,6,7,8,11,12,13,16,17,18,21,22,23)
colnames(fhsuadj2) <- c("y","ylo","yhi","x")



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p1c.pdf", height=6, width=4)
ggplot(fhsuadj2, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=c(rep(c("red","blue","green"),5))))+
geom_pointrange(size=1, shape=22)+
coord_flip() + 
geom_hline(aes(x=0), lty=2) + 
labs(x = "", y = "Beta for 1 SD increase in MG 18:2") + theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size = rel(1))) + guides(fill=FALSE) +
ylim(-0.6,+0.6) + 
ggtitle("MG 18:2") +
scale_x_discrete(breaks=c(2,7,12,17,22))
dev.off()


fhsuadj <- NULL
for (i in 1:nrow(resfhsU[[4]]))
{
	fhsuadj <- rbind(fhsuadj,rbind(resfhsU[[4]][i,],resfhsT[[4]][i,],resfhsP[[4]][i,]))
}

fhsuadj2 <- data.frame(as.numeric(fhsuadj[,2]),as.numeric(fhsuadj[,2])-1.96*as.numeric(fhsuadj[,3]),as.numeric(fhsuadj[,2])+1.96*as.numeric(fhsuadj[,3]))
fhsuadj2$x <- c(1,2,3,6,7,8,11,12,13,16,17,18,21,22,23)
colnames(fhsuadj2) <- c("y","ylo","yhi","x")



pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p1d.pdf", height=6, width=4)
ggplot(fhsuadj2, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=c(rep(c("red","blue","green"),5))))+
geom_pointrange(size=1, shape=22)+
coord_flip() + 
geom_hline(aes(x=0), lty=2) + 
labs(x = "", y = "Beta for 1 SD increase in SM 28:1") + theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size = rel(1))) + guides(fill=FALSE) +
ylim(-0.6,+0.6) + 
ggtitle("SM 28:1")+
scale_x_discrete(breaks=c(2,7,12,17,22))
dev.off()






#####################################################################################################################
## 2. FIGURE 2, PANEL B - ASSOCIATION BETWEEN 4 metabolites and inflammatory markers and subclinical CVD in PIVUS ###
#####################################################################################################################



### WITHOUT FHS RISK FACTORS ###

lysoPCL <- c("PE_cer","MG18_2","LPC18_2","LPC18_1")
## Association
RESuadj <- list(NA)
for (k in lysoPCL)
{
	temp <- NULL
	for (i in c("logcrp","icam1wlab","vcam1","lneselectin","lnil6","lntnfa","mcp1","ln.LP.PLA2",
	"oxldl","lnhomocyst","cd","gluratio",	"lvmim27","rwt","lnedv","lnfmd","lnccadistensibility","imtfmeansindx2","lnfibrinogen","lntpaag","lnpai1","lnddimer"
	))
	{
			st_r <- scale(as.numeric(mPIVUSF[,i]))
			mod <- lm(st_r~mPIVUSF$kvinna1+scale(as.numeric(mPIVUSF[,k])))
			temp <- rbind(temp,c(i,coefficients(mod)[3],coefficients(summary(mod))[3,2],coefficients(summary(mod))[3,4]))
	}
RESuadj[[which(k==lysoPCL)]] <- temp
}

### WITH FHS RISK FACTORS ###

lysoPCL <- c("PE_cer","MG18_2","LPC18_2","LPC18_1")
## Association
RESadj <- list(NA)
for (k in lysoPCL)
{
	temp <- NULL
	for (i in c("logcrp","icam1wlab","vcam1","lneselectin","lnil6","lntnfa","mcp1","ln.LP.PLA2",
	"oxldl","lnhomocyst","cd","gluratio",	"lvmim27","rwt","lnedv","lnfmd","lnccadistensibility","imtfmeansindx2","lnfibrinogen","lntpaag","lnpai1","lnddimer"
	))
	{
			st_r <- scale(as.numeric(mPIVUSF[,i]))
			mod <- lm(st_r~mPIVUSF$kvinna1+scale(as.numeric(mPIVUSF[,k]))+mPIVUSF$ldl+mPIVUSF$manuelltsbp+ mPIVUSF$hdl + mPIVUSF$bmi + mPIVUSF$smoke01 + mPIVUSF$diab+ mPIVUSF$logtg + mPIVUSF$antihyp)
			temp <- rbind(temp,c(i,coefficients(mod)[3],coefficients(summary(mod))[3,2],coefficients(summary(mod))[3,4]))
	}
RESadj[[which(k==lysoPCL)]] <- temp
}




y <- t(sapply(RESuadj, function(x) -log10(as.numeric(x[,4]))))
rownames(y) <- c("SM 28:1","MG 18:2","LysoPC 18:2","LysoPC 18:1")
colnames(y) <- c("Log C-reactive protein","ICAM-1","VCAM-1","Log E-selectin","Log Interleukin-6","Log Interferon alfa","MCP-1","Log Lp-PLA2 activity","Oxidized LDL","Log Homocystein","Conjugated dienes","Ox/Redox Glutathione ratio","Left ventricular mass index","Relative wall thickness", "Log End-diastolic volume", "Log Flow-mediated dilation", "Log CCA distensibility","Log Intima media thickness","Log Fibrinogen","Log Tissue plasminogen activator","Log Plasminogen activator inhibitor-1", "Log D-dimer")

ybeta <- t(sapply(RESuadj, function(x) as.numeric(x[,2])))

grid  <- expand.grid(x=seq(1:ncol(y)), y=seq(1:nrow(y)))
grid$z <- as.numeric(t(y))
grid$sign <- ifelse(as.numeric(t(ybeta))>0,"+","-")

library(lattice)
library(gplots)

pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p2a.pdf")
levelplot(z~x*y,grid, col.regions=rev(colorRampPalette(c("red","green", "aliceblue"), bias=0.5,interpolate="linear")(40)), main="-log10(P-value)",xlab="",ylab="",at=c(0,1.30,seq(2,15)), colorkey=list(space="top"),scales=list(x=list(rot=90,at=seq(1:ncol(y)),labels=colnames(y)),y = list(at=seq(1:nrow(y)),labels=rownames(y))),aspect="iso",
panel=function(...) {
                       arg <- list(...)
                       panel.levelplot(...)
                       panel.text(grid$x, grid$y,grid$sign)})
dev.off()





y <- t(sapply(RESadj, function(x) -log10(as.numeric(x[,4]))))
rownames(y) <- c("SM 28:1","MG 18:2","LysoPC 18:2","LysoPC 18:1")
colnames(y) <- c("Log C-reactive protein","ICAM-1","VCAM-1","Log E-selectin","Log Interleukin-6","Log Interferon alfa","MCP-1","Log Lp-PLA2 activity","Oxidized LDL","Log Homocystein","Conjugated dienes","Ox/Redox Glutathione ratio","Left ventricular mass index","Relative wall thickness", "Log End-diastolic volume", "Log Flow-mediated dilation", "Log CCA distensibility","Log Intima media thickness","Log Fibrinogen","Log Tissue plasminogen activator","Log Plasminogen activator inhibitor-1", "Log D-dimer")

ybeta <- t(sapply(RESadj, function(x) as.numeric(x[,2])))

grid  <- expand.grid(x=seq(1:ncol(y)), y=seq(1:nrow(y)))
grid$z <- as.numeric(t(y))
grid$sign <- ifelse(as.numeric(t(ybeta))>0,"+","-")

library(lattice)
library(gplots)

pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Fig_2_p2b.pdf")
levelplot(z~x*y,grid, col.regions=rev(colorRampPalette(c("red","green", "aliceblue"), bias=0.5,interpolate="linear")(40)), main="-log10(P-value)",xlab="",ylab="",at=c(0,1.30,3.32,seq(2,15)), colorkey=list(space="top"),scales=list(x=list(rot=90,at=seq(1:ncol(y)),labels=colnames(y)),y = list(at=seq(1:nrow(y)),labels=rownames(y))),aspect="iso",
panel=function(...) {
                       arg <- list(...)
                       panel.levelplot(...)
                       panel.text(grid$x, grid$y,grid$sign)})
dev.off()

