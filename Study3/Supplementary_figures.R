#----------------------------------------------
# Filename: Supplementary_figures.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 
# Purpose: Analysis used in the paper. Supplementary figures
# Note: 
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small), LPC18_1_metaanalysis1.tbl
# Data created: Suppl_figure1a.pdf Suppl_figure1b.pdf
#              for_locus_plot
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/

### FUCNTION ###

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



ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
 
  library(ggplot2)  
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
 
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
 
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
 
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
 
    pl <- ggplot(dat, aes(x = time, y = surv)) + 
      xlab(xlab) + ylab(ylab) + ggtitle(main) + 
      geom_step(col = col, lty = lty.est)
 
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
 
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                       col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations') 
    } else(pl)
 
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
 
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
 
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
 
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]), 
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1)) 
    }
 
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
 
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) + 
      xlab(xlab) + ylab(ylab) + ggtitle(main) + 
      geom_step(aes(col = group, lty = group))
 
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
 
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
 
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
 
    pl <- pl + line
 
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
               geom_step(aes(y = low,lty = group), col = surv.col)}   
    } else {pl}
 
 
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations') 
    } else(pl)
 
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl) 
    pl
  } 
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main) 
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}
###############
#### INDEX ####
### 1. SUPPLEMENTARY FIGURE 1 - Exploring interaction with age and linearity for MG and LysoPC
###                             (in the end we will keep only Suppl_figure1b.pdf and Suppl_figure1e.pdf)
### 2. SUPPLEMENTARY FIGURE 3 - LOCUS ZOOM PLOT



##################################################################################################
## 1. SUPPLEMENTARY FIGURE 1. Panel A - Exploring interaction with age and linearuty for LysoPC ##
##################################################################################################

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


#### LOAD ULSAM ####
library(survival)

load("/home/andrea/glob/alignment_ulsam_small/Results/Final_datasets/step4.Rdata")

write.csv(metabo_p[,c("pat","time","double_")], file="ulsam_key.csv", row.names=F)

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




##############################################
### AGE-VARYING EFFECT OF  LysoPC 18:2 #######
##############################################


library(survival)
library(pspline)
library(ggplot2)
library(grid)


##### LysoPC 18:2 ######


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi) & !is.na(metabo_sub_nopchdT$tg),]

metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

metabo_sub_nopchdT4$x <- scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671))


metabo_sub_nopchdT4 <- metabo_sub_nopchdT4[,c("x","sex","age","sbp","bmi","smoke01","antihyp","ldl","hdl","diab_baseline_new","surv_chd","incchd","weights","tg")]


mod <- coxph(Surv(surv_chd,incchd) ~ sex+bs(age)+x+bs(age):x, data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)

dfb<-as.matrix(resid(mod,type="dfbeta"))

# Recalculate the correct SE
strata_na <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

gamma<-matrix(0,dim(mod$var)[1],dim(mod$var)[2])
	for (s in 1:max(strata))
		{
				indst<-(1:length(metabo_sub_nopchdT4$surv_chd))[strata_na==s]
				m <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])
				n <- as.numeric(table(metabo_sub_nopchdT4$weights)[s+1])*as.numeric(names(table(metabo_sub_nopchdT4$weights)[s+1]))
				if (m>1) gamma<-gamma+(1-m/n)*m*var(dfb[indst,])
		} 
adjvar<-mod$var+gamma

m <- 300
newdata0 <- data.frame(
age=seq(min(metabo_sub_nopchdT4$age),max(metabo_sub_nopchdT4$age),length=m),
x=0,
sex=mean(metabo_sub_nopchdT4$sex))

newdata1 <- data.frame(
age=seq(min(metabo_sub_nopchdT4$age),max(metabo_sub_nopchdT4$age),length=m),
x=1,
sex=mean(metabo_sub_nopchdT4$sex))


X0 <- X.coxph(mod,newdata=newdata0)
X1 <- X.coxph(mod,newdata=newdata1)
lp <- predict(mod,newdata1) - predict(mod,newdata0)
se <- sqrt(diag((X1-X0) %*% adjvar %*% t(X1-X0)))


pl <- data.frame(MX=seq(min(metabo_sub_nopchdT4$age),max(metabo_sub_nopchdT4$age),length=m),MY=exp(lp),UY=exp(lp+1.96*se),LY=exp(lp-1.96*se))

library(grid)


p <- ggplot(pl, aes(MX)) + 
	ylim(0,1.5) + xlim(55,90) + 
	geom_hline(yintercept= 1 , colour="black" , linetype=2 ,lwd = 0.4)+
	geom_line(aes(y=MY), colour="blue") + 
  geom_ribbon(aes(ymin=UY, ymax=LY), alpha=0.2) +
	labs(x = "Age", y = "HR for 1 SD increase in LysoPC 18:2")
	
	for (i in metabo_sub_nopchdT4$age[metabo_sub_nopchdT4$incchd==1])
	{
		p <- p + annotation_custom(grob=linesGrob(),xmin=i,xmax=i,ymin=-0.2,ymax=0.05)
		print(i)
	}
	
pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Suppl_figure1a.pdf")
p
dev.off()





################################################################################
## 2. SUPPLEMENTARY FIGURE 2. Panel B - Survival curves for the 4 metabolites ##
################################################################################


metabo_sub_nopchdT4 <-  metabo_sub_nopchdT[ !is.na(metabo_sub_nopchdT$hdl) & !is.na(metabo_sub_nopchdT$ldl) & !is.na(metabo_sub_nopchdT$smoke01) & !is.na(metabo_sub_nopchdT$sbp) & !is.na(metabo_sub_nopchdT$antihyp) & !is.na(metabo_sub_nopchdT$diab_baseline_new) & !is.na(metabo_sub_nopchdT$bmi) & !is.na(metabo_sub_nopchdT$tg),]

metabo_sub_nopchdT4$stratum <- ifelse(metabo_sub_nopchdT4$incchd==1,0,metabo_sub_nopchdT4$stratum)

metabo_sub_nopchdT4$weights <- ifelse(metabo_sub_nopchdT4$stratum==0,1,
ifelse(metabo_sub_nopchdT4$stratum==1,metabo_sub_nopchdT4$strata1_n[1]/sum(metabo_sub_nopchdT4$stratum==1),
ifelse(metabo_sub_nopchdT4$stratum==2,metabo_sub_nopchdT4$strata2_n[1]/sum(metabo_sub_nopchdT4$stratum==2),
ifelse(metabo_sub_nopchdT4$stratum==3,metabo_sub_nopchdT4$strata3_n[1]/sum(metabo_sub_nopchdT4$stratum==3),metabo_sub_nopchdT4$strata4_n[1]/sum(metabo_sub_nopchdT4$stratum==4)))))

strata <- as.numeric(as.factor(metabo_sub_nopchdT4$weights))-1

## LysoPC 18:2 ###
vari1 <- scale(as.numeric(metabo_sub_nopchdT4$M520.340T346.671))
metabo_sub_nopchdT4$LPC182 <- cut(vari1,quantile(vari1,seq(0,1,length.out=4)))

## MG182 18:2 ###
vari2 <- scale(as.numeric(metabo_sub_nopchdT4$M393.231T384.578))
metabo_sub_nopchdT4$MG182 <- cut(vari2,quantile(vari1,seq(0,1,length.out=4)))

## LysoPC 18:1 ###
vari3 <- scale(as.numeric(metabo_sub_nopchdT4$M522.356T382.722))
metabo_sub_nopchdT4$LPC181 <- cut(vari3,quantile(vari1,seq(0,1,length.out=4)))

### PE-cer d34:1 ###
vari4 <- scale(as.numeric(metabo_sub_nopchdT4$M661.528T580.694))
metabo_sub_nopchdT4$PEcer  <- cut(vari4,quantile(vari4,seq(0,1,length.out=4)))


metabo_sub_nopchdT4 <- metabo_sub_nopchdT4[,c("sex","age","sbp","bmi","smoke01","antihyp","ldl","hdl","diab_baseline_new","surv_chd","incchd","weights","LPC182","MG182","LPC181","PEcer","tg")]


modLPC182 <- coxph(Surv(surv_chd/365.25,incchd) ~ sex+age+sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+LPC182+log(tg), data =metabo_sub_nopchdT4[metabo_sub_nopchdT4$age>70,],weights=metabo_sub_nopchdT4$weights[metabo_sub_nopchdT4$age>70])
modLMG182 <- coxph(Surv(surv_chd/365.25,incchd) ~ sex+age+sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+MG182+log(tg), data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)
modLPC181 <- coxph(Surv(surv_chd/365.25,incchd) ~ sex+age+sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+LPC181+log(tg), data =metabo_sub_nopchdT4[metabo_sub_nopchdT4$age>70,],weights=metabo_sub_nopchdT4$weights[metabo_sub_nopchdT4$age>70])
modPEcer <- coxph(Surv(surv_chd/365.25,incchd) ~ sex+age+sbp+bmi+smoke01+antihyp+ldl+hdl+diab_baseline_new+PEcer+log(tg), data =metabo_sub_nopchdT4,weights=metabo_sub_nopchdT4$weights)


ndLPC181 <- data.frame(
LPC181=levels(metabo_sub_nopchdT4$LPC181),
age=77,
sex=1,
sbp=150,
bmi=26,
smoke01=1,
antihyp=0,
ldl=2.6,
tg=1.7,
hdl=1.293,
diab_baseline_new=0)


ndLPC182 <- data.frame(
LPC182=levels(metabo_sub_nopchdT4$LPC182),
age=77,
sex=1,
sbp=150,
bmi=26,
smoke01=1,
antihyp=0,
ldl=2.6,
tg=1.7,
hdl=1.293,
diab_baseline_new=0)


ndMG182 <- data.frame(
MG182=levels(metabo_sub_nopchdT4$MG182),
age=77,
sex=1,
sbp=150,
bmi=26,
smoke01=1,
antihyp=0,
ldl=2.6,
tg=1.7,
hdl=1.293,
diab_baseline_new=0)


ndPEcer <- data.frame(
PEcer=levels(metabo_sub_nopchdT4$PEcer),
age=77,
sex=1,
sbp=150,
bmi=26,
smoke01=1,
antihyp=0,
ldl=2.6,
tg=1.7,
hdl=1.293,
diab_baseline_new=0)



forsurvplotLPC181 <- survfit(modLPC181,ndLPC181)
forsurvplotLPC182 <- survfit(modLPC182,ndLPC182)
forsurvplotMG182 <- survfit(modLMG182,ndMG182)
forsurvplotPEcer <- survfit(modPEcer,ndPEcer)

pdf("/home/andrea/glob/alignment_pivus_small/Results/chd/paper/Suppl_figure1b.pdf")
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(forsurvplotLPC181, ylim=c(0.80,1), mark.time=T, col=c("red","yellow","green"),mark=1, xlab="Years", ylab="Survival",main="LysoPC 18:1")
legend("bottomleft",inset=0.05,c("1st tertile","2nd tertile","3rd tertile"), col=c("red","yellow","green"), pch=16,bty="n")
plot(forsurvplotLPC182, ylim=c(0.80,1), mark.time=T, col=c("red","yellow","green"),mark=1, xlab="Years", ylab="Survival",main="LysoPC 18:2")
plot(forsurvplotMG182, ylim=c(0.80,1), mark.time=T, col=c("red","yellow","green"),mark=1, xlab="Years", ylab="Survival",main="MG 18:2")
plot(forsurvplotPEcer, ylim=c(0.80,1), mark.time=T, col=c("red","yellow","green"),mark=1, xlab="Years", ylab="Survival",main="SM 28:1")
dev.off()



#######################################################
##### 2. SUPPLEMENTARY FIGURE 3 - LOCUS ZOOM PLOT #####
#######################################################


### Save data for locus-zoom

ph <- "LPC18_1"
cc <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
t <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)


### READ SNP ANNOTATION FILE ###
write.table(t[t[,4] >= 0.03 & t[,4] <= 0.97,c(1,10)], file="/proj/b2011036/nobackup/for_locus_plot", row.names=F, quote=F)

### The run LOCUS-ZOOM ###





