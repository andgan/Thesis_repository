#----------------------------------------------
# Filename: Figures.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 06SEP2014
# Updated: 
# Purpose: Program to create the figures reported in the paper
# Note: 
#-----------------------------------------------
# Data used: univariateM.xlsx univariateF.xlsx loeEF.Rdata hlexpEF.Rdata loeEM.Rdata hlexpEM.Rdata DISCCALM.Rdata MM095.Rdata FF095.Rdata All_riskF.Rdata All_riskM.Rdata AMV1EX.Rdata AFV1EX.Rdata male_2009_11.csv female_2009_11.csv AM1 AF1
# Data created: Figure1 to Suppl Figure 13.pdf
#-----------------------------------------------
# OP: R 3.1.0
#-----------------------------------------------*/

## Load function ##
library(gdata)
library(ggplot2)
library(plyr)
library(survival)
library(scales)

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
# READ MANUALLY PROCESSED FILE (manually procssed version from univariateM.tab) #
resM <- read.xls("/proj/b2011036/uk.biobank/univariateM.xlsx", header=T, stringsAsFactor=F)
resF <- read.xls("/proj/b2011036/uk.biobank/univariateF.xlsx", header=T, stringsAsFactor=F)

cbbPalette <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","BROWN")

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
		
		K <- rbind(K,c(f,unique(t$R2.for.association.with.age),unique(t$C4),unique(t$Measurement.class),k))
	}	
	
	assign(paste0("MCIN",dis),K)	
}		
		
		

### PROCESS DATA IN WOMEN ###
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
		
		K <- rbind(K,c(f,unique(t$R2.for.association.with.age),unique(t$C4),unique(t$Measurement.class),k))
	}	
	
	assign(paste0("FCIN",dis),K)	
}		
		


################
### FIGURE 1 ###
################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINall[!is.na(MCINall[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4 
resMcl2$Category42[-ind] <- ""
resMcl2$shape <- 19
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(-0.02,-0.02,1,1,1,1,-0.02,1,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0,0,0,0,0,0,0,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,70,0,0,30,30,70,0,30,30)

resMcl2$Code[order(-resMcl2$C.index)[1:10]]
resMcl2$Category42[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Figure1M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at,shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) + 
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()


# Female #
resFcl2 <- data.frame(FCINall[!is.na(FCINall[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index) 

ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape <- 19
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,1,-0.02,1,1,1,1,1,1,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,1.2,0,0,0,0,0,0,0,0)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Figure1F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at,shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) + 
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



################
### FIGURE 2 ###
################


## PREPARE DATA ##
# Male #
resMFcl <- merge(data.frame(MCINall[!is.na(MCINall[,5]),], stringsAsFactors = F),data.frame(FCINall[!is.na(FCINall[,5]),], stringsAsFactors = F), by="X1")

cbbPalette2 <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","BROWN")


resMFcl2 <- resMFcl[,c(1,2,3,4,5,9)] 
colnames(resMFcl2) <- c("Code", "R.2","Category4","Category2_at","C.indexM","C.indexF")
resMFcl2$R.2 <- as.numeric(resMFcl2$R.2) 
resMFcl2$C.indexM <- as.numeric(resMFcl2$C.indexM) 
resMFcl2$C.indexF <- as.numeric(resMFcl2$C.indexF) 

res <- abs(residuals(lm(resMFcl2$C.indexF ~ resMFcl2$C.indexM)))
ind <- order(-res)[1:10] 
resMFcl2$Category42 <- resMFcl2$Category4
resMFcl2$Category42[-ind] <- ""
resMFcl2$shape <- 19
resMFcl2$shape[ind] <- 17
resMFcl2$vjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$vjust[ind] <- c(0,-0.8,0,0.3,0,-0.5,0,0,0,1.3)
resMFcl2$hjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$hjust[ind] <- c(-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,1,-0.02,-0.02)
#resMFcl2$angle <- rep(0,1,nrow(resMFcl2))
#resMFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resMFcl2$Category42[ind]
resMFcl2$Code[ind]
resMFcl2$C.indexM[ind]
resMFcl2$C.indexF[ind]

pdf("/proj/b2011036/uk.biobank/figures/Figure2MF.pdf", width=10)
ggplot(resMFcl2,aes(x= C.indexM, y= C.indexF)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() +
scale_colour_manual(values=cbbPalette2,name = "Category") +
xlab("C-index including age - Men") + ylab("C-index including age -  Women") + stat_smooth(method="lm", se=FALSE) + scale_x_continuous(breaks = seq(0.67, 0.75, by = 0.02)) +  geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=0),size=3) + ylim(0.65, 0.73)
dev.off()




##############
## FIGURE 3 ## 
##############

# Create a healthy c-index dataset for merging #

## PREPARE DATA ##
# Male #


resMcl2 <- data.frame(MCINHealthy[!is.na(MCINHealthy[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

ind <- order(-resMcl2$C.index)[14:22]
resMcl2$Category42 <- resMcl2$Category4 
resMcl2$Category42[-ind] <- ""
resMcl2$shape <- 19
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(-0.02,1,1,1,1,-0.02,1,1,1,-0.02)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0,0,0,0,0,0,0,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,0,0,30,0,0,0,0)

resMcl2$Code[order(-resMcl2$C.index)[14:22]]
resMcl2$Category42[order(-resMcl2$C.index)[14:22]]


ann_text <- data.frame(R.2 = 0.0006,C.index = 0.702,lab = "SMOKING-RELATED QUESTIONS")


pdf("/proj/b2011036/uk.biobank/figures/Figure3M.pdf", width=10)
ggplot(resMcl2,aes(x=R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at,shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) + 
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x))) +
geom_text(data = ann_text,label = "Smoking-related questions", size=3) + ylim(0.677, 0.703)
dev.off()



# Female

resFcl2 <- data.frame(FCINHealthy[!is.na(FCINHealthy[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index) 

ind <- order(-resFcl2$C.index)[15:24]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape <- 19
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,1,1,1,1,1,1,1,1,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,0,0,0.4,0,0,0,0,0,0)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,30,30,30,30,0,30,30)

resFcl2$Category42[ind]
resFcl2$Code[ind]

ann_text <- data.frame(R.2 = 0.0004,C.index = 0.682,lab = "SMOKING-RELATED QUESTIONS")


pdf("/proj/b2011036/uk.biobank/figures/Figure3F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at,shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) + 
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x))) +
geom_text(data = ann_text,label = "Smoking-related questions", size=3) + ylim(0.660, 0.685)
dev.off()




#######################################
### SUPPLEMENTARY FIGURE 1 - CANCER ###
#######################################


## PREPARE DATA ##
# Male #

resMcl2 <- data.frame(MCINCA[!is.na(MCINCA[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(1,1,-0.02,-0.02,1,1,-0.02,1,-0.02,+1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0.0,0,0,0,0,0,0.8,0.0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure1M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()




# Female #
resFcl2 <- data.frame(FCINCA[!is.na(FCINCA[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)

resFcl2$shape <- 19
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,1,1,1,1,1,1,1,1,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,1.2,0,0,-0.7,1.2,0,-0.7,0,0)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure1F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



###############################################
### SUPPLEMENTARY FIGURE 2 - CARDIOVASCULAR ###
###############################################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINCVD[!is.na(MCINCVD[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(1,-0.02,-0.02,1,1,1,1,1,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0.0,0,0,0.2,0,0,0.0,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,0,0,0,0,30,30,0)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure2M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



# Female #
resFcl2 <- data.frame(FCINCVD[!is.na(FCINCVD[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)

resFcl2$shape<- 19
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(-0.02,1,1,1,-0.02,-0.02,+1,-0.02,1,-0.02)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,0,0.5,0,0,0,-0.5,0,0.5,1.5)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,50,50,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure2F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



############################################
### SUPPLEMENTARY FIGURE 3 - RESPIRATORY ###
############################################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINRE[!is.na(MCINRE[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(-0.02,-0.02,1,-0.02,1,1,1,1,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0.0,0,0,-0.5,0.2,0,0.0,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,30,30)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure3M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



# Female #
resFcl2 <- data.frame(FCINRE[!is.na(FCINRE[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)

resFcl2$shape <- 19
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,-0.02,-0.02,-0.02,1,1,1,-0.02,-0.02,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,0,0,0,-0.8,0,0,0,0,0)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure3F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x))) + ylim(0.66,0.83)
dev.off()





##########################################
### SUPPLEMENTARY FIGURE 4 - DIGESTIVE ###
##########################################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINDG[!is.na(MCINDG[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(1,1,1,1,1,-0.02,1,1,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0.5,0.0,0,0,0.5,0,0.5,0,1)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,30,0,30,30,30,0,30)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure4M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



# Female #
resFcl2 <- data.frame(FCINDG[!is.na(FCINDG[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)

resFcl2$shape <- 19
resFcl2$indf2 <- factor(resFcl2$indf, labels = c("All variables", "Categories"))
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(-0.02,1,1,1,1,1,1,1,1,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,0,0,0,-0.8,0,0.5,0,0,0.8)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure4F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



################################################
### SUPPLEMENTARY FIGURE 5 - EXTERNAL CAUSES ###
################################################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINEC[!is.na(MCINEC[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(1,1,1,1,1,1,1,1,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0.5,0.8,0,0,0,0,0,0.3,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,30,30,30,30,30,30,30)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure5M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()




# Female #
resFcl2 <- data.frame(FCINEC[!is.na(FCINEC[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)

resFcl2$shape <- 19
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,1,1,1,1,1,1,1,1,1)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(0,0,0,0,-0.5,0,1.2,0,0,1)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,0,0,0,0,30,0,0,30,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure5F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()





##############################################
### SUPPLEMENTARY FIGURE 6 - OTHER DISEASE ###
##############################################


## PREPARE DATA ##
# Male #
resMcl2 <- data.frame(MCINOD[!is.na(MCINOD[,5]),], stringsAsFactors = F)
colnames(resMcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resMcl2$R.2 <- as.numeric(resMcl2$R.2) 
resMcl2$C.index <- as.numeric(resMcl2$C.index) 

resMcl2$shape <- 19
ind <- order(-resMcl2$C.index)[1:10]
resMcl2$Category42 <- resMcl2$Category4
resMcl2$Category42[-ind] <- ""
resMcl2$shape[ind] <- 17
resMcl2$hjust <- rep(0,1,nrow(resMcl2))
resMcl2$hjust[ind] <- c(1,-0.02,-0.02,-0.02,1,1,-0.02,-0.02,1,1)
resMcl2$vjust <- rep(0,1,nrow(resMcl2))
resMcl2$vjust[ind] <- c(0,0,0.0,0.3,0,0,0.8,0,0.8,0)
resMcl2$angle <- rep(0,1,nrow(resMcl2))
resMcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)


resMcl2$Category42[order(-resMcl2$C.index)[1:10]]
resMcl2$Code[order(-resMcl2$C.index)[1:10]]


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure6M.pdf", width=10)
ggplot(resMcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



# Female #
resFcl2 <- data.frame(FCINOD[!is.na(FCINOD[,5]),], stringsAsFactors = F)
colnames(resFcl2) <- c("Code", "R.2","Category4","Category2_at","C.index")
resFcl2$R.2 <- as.numeric(resFcl2$R.2) 
resFcl2$C.index <- as.numeric(resFcl2$C.index)


resFcl2[resFcl2$R.2==min(resFcl2$R.2),]

resFcl2$shape <- 19
ind <- order(-resFcl2$C.index)[1:10]
resFcl2$Category42 <- resFcl2$Category4
resFcl2$Category42[-ind] <- ""
resFcl2$shape[ind] <- 17
resFcl2$hjust <- rep(0,1,nrow(resFcl2))
resFcl2$hjust[ind] <- c(1,1,1,-0.02,1,-0.02,-0.02,1,1,-0.02)
resFcl2$vjust <- rep(0,1,nrow(resFcl2))
resFcl2$vjust[ind] <- c(-0.8,0,0,0,0,0.8,0,0,0,0)
resFcl2$angle <- rep(0,1,nrow(resFcl2))
resFcl2$angle[ind] <- c(0,5,15,0,0,0,0,0,0,0)

resFcl2$Category42[order(-resFcl2$C.index)[1:10]]
resFcl2$Code[order(-resFcl2$C.index)[1:10]]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure6F.pdf", width=10)
ggplot(resFcl2,aes(x= R.2, y= C.index)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() + 
scale_colour_manual(values=cbbPalette,name = "Category") +
xlab(expression("Association with age (" ~ italic(R)^2 ~ ")")) +
ylab("C-index including age") + 
geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + scale_x_continuous(trans="log10", breaks = c(10^-8,10^-6, 10^-4,10^-2,10^-1), labels = trans_format("log10", math_format(10^.x)))
dev.off()



########################################
### SUPPLEMENTARY FIGURE 7 - CANCER  ###
########################################


## PREPARE DATA ##
resMFcl <- merge(data.frame(MCINCA[!is.na(MCINCA[,5]),], stringsAsFactors = F),data.frame(FCINCA[!is.na(FCINCA[,5]),], stringsAsFactors = F), by="X1")

cbbPalette2 <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","BROWN")


resMFcl2 <- resMFcl[,c(1,2,3,4,5,9)] 
colnames(resMFcl2) <- c("Code", "R.2","Category4","Category2_at","C.indexM","C.indexF")
resMFcl2$R.2 <- as.numeric(resMFcl2$R.2) 
resMFcl2$C.indexM <- as.numeric(resMFcl2$C.indexM) 
resMFcl2$C.indexF <- as.numeric(resMFcl2$C.indexF) 

resMFcl2$shape <- 19
res <- abs(residuals(lm(resMFcl2$C.indexF ~ resMFcl2$C.indexM)))
ind <- order(-res)[1:10] 
resMFcl2$Category42 <- resMFcl2$Category4
resMFcl2$Category42[-ind] <- ""
resMFcl2$shape[ind] <- 17
resMFcl2$vjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$vjust[ind] <- c(0,-0.4,0,0.8,0,0.2,0,0.8,0,0.4)
resMFcl2$hjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$hjust[ind] <- c(1,1,0,0,0,0,0,0,0,0)
resMFcl2$angle <- rep(0,1,nrow(resMFcl2))
resMFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,30,0)

resMFcl2$Category42[ind]
resMFcl2$Code[ind]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure7MF.pdf", width=10)
ggplot(resMFcl2,aes(x= C.indexM, y= C.indexF)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() +
scale_colour_manual(values=cbbPalette2,name = "Category") +
xlab("C-index  including age - Men") + ylab("C-index  including age - Women") + stat_smooth(method="lm", se=FALSE) + scale_x_continuous(breaks = seq(0.69, 0.77, by = 0.02)) +  geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + ylim(0.64, 0.76)
dev.off()





###############################################
### SUPPLEMENTARY FIGURE 8 - CARDIOVASCULAR ###
###############################################


## PREPARE DATA ##
resMFcl <- merge(data.frame(MCINCVD[!is.na(MCINCVD[,5]),], stringsAsFactors = F),data.frame(FCINCVD[!is.na(FCINCVD[,5]),], stringsAsFactors = F), by="X1")

cbbPalette2 <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","BROWN")


resMFcl2 <- resMFcl[,c(1,2,3,4,5,9)] 
colnames(resMFcl2) <- c("Code", "R.2","Category4","Category2_at","C.indexM","C.indexF")
resMFcl2$R.2 <- as.numeric(resMFcl2$R.2) 
resMFcl2$C.indexM <- as.numeric(resMFcl2$C.indexM) 
resMFcl2$C.indexF <- as.numeric(resMFcl2$C.indexF) 


resMFcl2$shape <- 19
res <- abs(residuals(lm(resMFcl2$C.indexF ~ resMFcl2$C.indexM)))
ind <- order(-res)[14:23] 
resMFcl2$Category42 <- resMFcl2$Category4
resMFcl2$Category42[-ind] <- ""
resMFcl2$shape[ind] <- 17
resMFcl2$vjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$vjust[ind] <- c(0,0,0,0,0,0.7,0,0,0,0)
resMFcl2$hjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$hjust[ind] <- c(0,0,0,0,0,0,0,0,0,0)
resMFcl2$angle <- rep(0,1,nrow(resMFcl2))
resMFcl2$angle[ind] <- c(0,30,30,0,0,0,50,0,0,30)

resMFcl2$Category42[ind]
resMFcl2$Code[ind]

ann_text <- data.frame(C.indexM = 0.695,C.indexF = 0.725,lab = "SMOKING-RELATED QUESTIONS")


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure8MF.pdf", width=10)
ggplot(resMFcl2,aes(x= C.indexM, y= C.indexF)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() +
scale_colour_manual(values=cbbPalette2,name = "Category") +
xlab("C-index  including age - Men") + ylab("C-index  including age - Women") + stat_smooth(method="lm", se=FALSE) + scale_x_continuous(breaks = seq(0.68, 0.77, by = 0.02)) +  geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + ylim(0.68, 0.75) + geom_text(data = ann_text,label = "Smoking-related questions", size=3)
dev.off()





############################################
### SUPPLEMENTARY FIGURE 9 - RESPIRATORY ###
############################################



resMFcl <- merge(data.frame(MCINRE[!is.na(MCINRE[,5]),], stringsAsFactors = F),data.frame(FCINRE[!is.na(FCINRE[,5]),], stringsAsFactors = F), by="X1")

cbbPalette2 <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","BROWN")


resMFcl2 <- resMFcl[,c(1,2,3,4,5,9)] 
colnames(resMFcl2) <- c("Code", "R.2","Category4","Category2_at","C.indexM","C.indexF")
resMFcl2$R.2 <- as.numeric(resMFcl2$R.2) 
resMFcl2$C.indexM <- as.numeric(resMFcl2$C.indexM) 
resMFcl2$C.indexF <- as.numeric(resMFcl2$C.indexF) 

resMFcl2$shape <- 19
res <- abs(residuals(lm(resMFcl2$C.indexF ~ resMFcl2$C.indexM)))
ind <- order(-res)[1:10] 
resMFcl2$Category42 <- resMFcl2$Category4
resMFcl2$Category42[-ind] <- ""
resMFcl2$shape[ind] <- 17
resMFcl2$vjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$vjust[ind] <- c(0,0,0,0,0,0.5,1,-0.6,0,0)
resMFcl2$hjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$hjust[ind] <- c(0,0,0,0,1.05,0,0,0,0,1)
resMFcl2$angle <- rep(0,1,nrow(resMFcl2))
resMFcl2$angle[ind] <- c(0,0,0,0,30,10,0,30,30,30)

resMFcl2$Category42[ind]
resMFcl2$Code[ind]



pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure9MF.pdf", width=10)
ggplot(resMFcl2,aes(x= C.indexM, y= C.indexF)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() +
scale_colour_manual(values=cbbPalette2,name = "Category") +
xlab("C-index  including age - Men") + ylab("C-index  including age - Women") + stat_smooth(method="lm", se=FALSE)  + scale_x_continuous(breaks = seq(0.72, 0.85, by = 0.02)) +  geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + ylim(0.67, 0.83)
dev.off()






###########################################
### SUPPLEMENTARY FIGURE 10 - DIGESTIVE ###
###########################################


resMFcl <- merge(data.frame(MCINDG[!is.na(MCINDG[,5]),], stringsAsFactors = F),data.frame(FCINDG[!is.na(FCINDG[,5]),], stringsAsFactors = F), by="X1")

cbbPalette2 <- c("blue","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","BROWN")


resMFcl2 <- resMFcl[,c(1,2,3,4,5,9)] 
colnames(resMFcl2) <- c("Code", "R.2","Category4","Category2_at","C.indexM","C.indexF")
resMFcl2$R.2 <- as.numeric(resMFcl2$R.2) 
resMFcl2$C.indexM <- as.numeric(resMFcl2$C.indexM) 
resMFcl2$C.indexF <- as.numeric(resMFcl2$C.indexF) 

resMFcl2$shape <- 19
res <- abs(residuals(lm(resMFcl2$C.indexF ~ resMFcl2$C.indexM)))
ind <- order(-res)[1:10] 
resMFcl2$Category42 <- resMFcl2$Category4
resMFcl2$Category42[-ind] <- ""
resMFcl2$shape[ind] <- 17
resMFcl2$vjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$vjust[ind] <- c(0,0,0,0,0,0.8,0,0,0,-0.7)
resMFcl2$hjust <- rep(0,1,nrow(resMFcl2))
resMFcl2$hjust[ind] <- c(1,1,1,0,0,0,0,1,0,1)
resMFcl2$angle <- rep(0,1,nrow(resMFcl2))
resMFcl2$angle[ind] <- c(0,0,0,0,0,0,0,0,0,0)

resMFcl2$Category42[ind]
resMFcl2$Code[ind]
resMFcl2$C.indexM[ind]
resMFcl2$C.indexF[ind]

pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure10MF.pdf", width=10)
ggplot(resMFcl2,aes(x= C.indexM, y= C.indexF)) + 
geom_point(aes(colour=Category2_at, shape=shape), size=2) + scale_shape_identity() +
scale_colour_manual(values=cbbPalette2,name = "Category") +
xlab("C-index  including age - Men") + ylab("C-index  including age - Women") + stat_smooth(method="lm", se=FALSE) + scale_x_continuous(breaks = seq(0.55, 0.77, by = 0.02)) +  geom_text(aes(label=Category42,hjust=hjust, vjust=vjust, angle=angle),size=3) + ylim(0.63, 0.81)
dev.off()




######################################################
######  SUPPLEMENTARY FIGURE 11 - CALIBRATION  #######
######################################################


load("/proj/b2011036/uk.biobank/prediction_results/loeEF.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/hlexpEF.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/loeEM.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/hlexpEM.Rdata")


d <-  data.frame(x=c(loeEF[,1],loeEM[,1]), y=c(loeEF[,2],loeEM[,2]), cl=as.factor(c(rep(17,nrow(loeEF)), rep(16,nrow(loeEM)))))
d2 <-  data.frame(y=c(hlexpEF[,1],hlexpEM[,1]), x=c(hlexpEF[,2],hlexpEM[,2]), cl=as.factor(c(rep(17,nrow(hlexpEF)), rep(16,nrow(hlexpEM)))))


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure11.pdf", width=8)
ggplot() + 
geom_point(data=d2,aes(x=x, y=y,color=cl, group=cl,shape=cl), size=3) +
xlab("Predicted risk") + ylab("Observed risk") +
xlim(0,0.25) + ylim(0,0.25) +
geom_abline(intercept = 0, slope = 1) +
geom_line(data=d,aes(x=x, y=y,colour=cl, group=cl),linetype="dotted", size=1) + scale_colour_manual(values=c("blue","red"), labels = c("Females","Males"), name="") + scale_shape_manual(values=c(17,16), labels = c("Females","Males"), name="") 
dev.off()


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure11zoom.pdf", width=8)
ggplot() + 
geom_point(data=d2,aes(x=x, y=y,color=cl, group=cl,shape=cl), size=3) +
xlab("Predicted risk") + ylab("Observed risk") +
xlim(0,0.05) + ylim(0,0.05) +
geom_abline(intercept = 0, slope = 1) +
geom_line(data=d,aes(x=x, y=y,colour=cl, group=cl),linetype="dotted", size=1) + scale_colour_manual(values=c("blue","red"), labels = c("Females","Males"), name="") + scale_shape_manual(values=c(17,16), labels = c("Females","Males"), name="") 
dev.off()



##########################################################################
###### SUPPLEMENTARY FIGURE 12 - SURVIVAL IN DIFFERENT POPULATIONS #######
##########################################################################

load("/proj/b2011036/uk.biobank/imputation_results/AM1.Rdata")
load("/proj/b2011036/uk.biobank/imputation_results/AF1.Rdata")

set.seed(123)
idsplitM2 <- which(AM1$center%in%c("Glasgow","Edinburgh"))

AMT1 <- AM1[-unique(c(idsplitM2)),]
AMV1 <- AM1[idsplitM2,]

set.seed(123)
idsplitF2 <- which(AF1$center%in%c("Glasgow","Edinburgh"))

AFT1 <- AF1[-unique(c(idsplitF2)),]
AFV1 <- AF1[idsplitF2,]

maxsurv <- 5

# Baseline survival at 5-years for a fitiscious population with same UK mortality but UK biobank age-distribution #

muktM <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/male_2009_11.csv")
muktF <- read.csv("/proj/b2011036/uk.biobank/uk_national_stat/female_2009_11.csv")

ALIVEM <- NULL
DEATHM <- NULL
ALIVEF <- NULL
DEATHF <- NULL
for (i in 40:70)
{
	# Age-specific 5 year survival
	M <- muktM[muktM[,1]%in%(i+maxsurv),4]/muktM[muktM[,1]%in%i,4]
	F <- muktF[muktF[,1]%in%(i+maxsurv),4]/muktF[muktF[,1]%in%i,4]
	
	aliveM <- as.numeric(M)*sum(AM1$age==i)
	deathM <- sum(AM1$age==i)-aliveM
	aliveF <- as.numeric(F)*sum(AF1$age==i)
	deathF <- sum(AF1$age==i)-aliveF

	# Number of alive and deaths for each age using the UK mortality
	ALIVEM <- c(ALIVEM,aliveM)
	DEATHM <- c(DEATHM,deathM)
	ALIVEF <- c(ALIVEF,aliveF)
	DEATHF <- c(DEATHF,deathF)
}


# Get the weibull distribution from UK biobank
KM <- survreg(Surv(surv,out)~1, data=AM1, dist="weibull",maxiter=100)
KF <- survreg(Surv(surv,out)~1, data=AF1, dist="weibull",maxiter=100)

# The deaths distribution is samples from the weibull distribution as in the UK biobank
rrM <- rweibull(10000000, shape=1/KM$scale, scale = exp(coef(KM)[1]))
rsM <- sample(rrM[rrM<5],sum(DEATHM))

rrF <- rweibull(10000000, shape=1/KF$scale, scale = exp(coef(KF)[1]))
rsF <- sample(rrF[rrF<5],sum(DEATHF))

# Set those that are alive at the end of follow-up
survoM <- c(rep(5,sum(ALIVEM)),rsM)
outoM <- c(rep(0,sum(ALIVEM)),rep(1,sum(DEATHM)))
modoM <- survfit(Surv(survoM,outoM)~1)

survoF <- c(rep(5,sum(ALIVEF)),rsF)
outoF <- c(rep(0,sum(ALIVEF)),rep(1,sum(DEATHF)))
modoF <- survfit(Surv(survoF,outoF)~1)


d <- data.frame(y=c(survfit(Surv(surv, out)~1,data=AMT1,se.fit=T)$surv,survfit(Surv(surv, out)~1,data=AMV1,se.fit=T)$surv,modoM$surv),
x=c(survfit(Surv(surv, out)~1,data=AMT1,se.fit=T)$time,survfit(Surv(surv, out)~1,data=AMV1,se.fit=T)$time,modoM$time),
U=c(survfit(Surv(surv, out)~1,data=AMT1,se.fit=T)$upper,survfit(Surv(surv, out)~1,data=AMV1,se.fit=T)$upper,modoM$surv), L=c(survfit(Surv(surv, out)~1,data=AMT1,se.fit=T)$lower,survfit(Surv(surv, out)~1,data=AMV1,se.fit=T)$lower,modoM$surv),
cl=c(rep(1,length(survfit(Surv(surv, out)~1,data=AMT1,se.fit=T)$surv)),rep(2,length(survfit(Surv(surv, out)~1,data=AMV1,se.fit=T)$surv)),rep(3,length(modoM$surv))))


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure12M.pdf", width=10)
ggplot(d) + 
geom_line(aes(x=x, y=y, group=cl, linetype=as.factor(cl)), color="black") + 
geom_ribbon(aes(ymin=L, ymax=U, x=x,group=cl), alpha=0.2) +
xlab("Time-in-study (in years)") + ylab("Survival") +
xlim(0,5) + ylim(0.95,1) + scale_linetype_manual(values=c(1,2,3), labels=c("Uk Biobank","Uk Biobank: \n Glasgow & Edinburgh centers", "Uk population"), name="")
dev.off()



d <- data.frame(y=c(survfit(Surv(surv, out)~1,data=AFT1,se.fit=T)$surv,survfit(Surv(surv, out)~1,data=AFV1,se.fit=T)$surv,modoF$surv),
x=c(survfit(Surv(surv, out)~1,data=AFT1,se.fit=T)$time,survfit(Surv(surv, out)~1,data=AFV1,se.fit=T)$time,modoF$time),
U=c(survfit(Surv(surv, out)~1,data=AFT1,se.fit=T)$upper,survfit(Surv(surv, out)~1,data=AFV1,se.fit=T)$upper,modoF$surv), L=c(survfit(Surv(surv, out)~1,data=AFT1,se.fit=T)$lower,survfit(Surv(surv, out)~1,data=AFV1,se.fit=T)$lower,modoF$surv),
cl=c(rep(1,length(survfit(Surv(surv, out)~1,data=AFT1,se.fit=T)$surv)),rep(2,length(survfit(Surv(surv, out)~1,data=AFV1,se.fit=T)$surv)),rep(3,length(modoF$surv))))


pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure12F.pdf", width=10)
ggplot(d) + 
geom_line(aes(x=x, y=y, group=cl, linetype=as.factor(cl)), color="black") + 
geom_ribbon(aes(ymin=L, ymax=U, x=x,group=cl), alpha=0.2) +
xlab("Time-in-study (in years)") + ylab("Survival") +
xlim(0,5) + ylim(0.95,1) + scale_linetype_manual(values=c(1,2,3), labels=c("Uk Biobank","Uk Biobank: \n Glasgow & Edinburgh centers", "Uk population"), name="")
dev.off()


################################################
### SUPPLEMENTARY FIGURE 13 - C-INDEX BY AGE ###
################################################

### BIOLOGICAL AGE ###

load("/proj/b2011036/uk.biobank/prediction_results/AMV1EX.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/All_riskM.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/AFV1EX.Rdata")
load("/proj/b2011036/uk.biobank/prediction_results/All_riskF.Rdata")


CC1AGM <- NULL
for (i in 40:70)
{
	cc1ag <- survConcordance(Surv(AMV1EX$surv[AMV1EX$age==i],AMV1EX$out[AMV1EX$age==i]==1)~All_riskM[AMV1EX$age==i])$concordance
	CC1AGM <- c(CC1AGM,cc1ag)
}

CC1AGF<- NULL
for (i in 40:70)
{
	cc1ag <- survConcordance(Surv(AFV1EX$surv[AFV1EX$age==i],AFV1EX$out[AFV1EX$age==i]==1)~All_riskF[AFV1EX$age==i])$concordance
	CC1AGF <- c(CC1AGF,cc1ag)
}

library("scales")

d <- data.frame(C=c(CC1AGM,CC1AGF),age=c(40:70,40:70),cl=c(rep(1,length(CC1AGM)),rep(2,length(CC1AGF))))

pdf("/proj/b2011036/uk.biobank/figures/Suppl_Figure13.pdf", width=10)
ggplot(d) + 
stat_smooth(aes(x=age, y=C, group=cl, fill=as.factor(cl),colour=as.factor(cl)), size=2, alpha = 0.3,span = 0.9) +
xlab("Age") + ylab("C-index") +
xlim(40,70) + scale_color_manual(values=c("blue","orange"), labels=c("Men","Women"), name="") + scale_fill_manual(values=c("blue","orange")) + guides(fill=FALSE) +  scale_y_continuous(limit=c(0.5,1),oob=squish)
dev.off()


## TO ADD IN THE PAPER ##

survConcordance(Surv(AFV1EX$surv[AFV1EX$age>40 & AFV1EX$age<50],AFV1EX$out[AFV1EX$age>40 & AFV1EX$age<50]==1)~All_riskF[AFV1EX$age>40 & AFV1EX$age<50])$concordance
survConcordance(Surv(AFV1EX$surv[AFV1EX$age>60 & AFV1EX$age<70],AFV1EX$out[AFV1EX$age>60 & AFV1EX$age<70]==1)~All_riskF[AFV1EX$age>60 & AFV1EX$age<70])$concordance


survConcordance(Surv(AMV1EX$surv[AMV1EX$age>40 & AMV1EX$age<50],AMV1EX$out[AMV1EX$age>40 & AMV1EX$age<50]==1)~All_riskM[AMV1EX$age>40 & AMV1EX$age<50])$concordance
survConcordance(Surv(AMV1EX$surv[AMV1EX$age>60 & AMV1EX$age<70],AMV1EX$out[AMV1EX$age>60 & AMV1EX$age<70]==1)~All_riskM[AMV1EX$age>60 & AMV1EX$age<70])$concordance
