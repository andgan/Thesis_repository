#----------------------------------------------
# Filename: Table2.R
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 12JAN2014
# Updated: 
# Purpose: Analysis used in the paper. Table 2
# Note: GWAS for 4 top-metabolites
#-----------------------------------------------
# Data used: step4.Rdata (from Twingene Small, Ulsam small, Pivus small)
# Data created: {ph}.metaanalysis_e_6.csv
#-----------------------------------------------
# OP: R 2.13.1, 
#-----------------------------------------------*/

########## 
### FUNCTIONS 
#########



## qq-plot
gg.qq = function(pvector, title=NULL, spartan=F) {
	library(ggplot2)
	o = -log10(sort(pvector,decreasing=F))
	#e = -log10( 1:length(o)/length(o) )
	e = -log10( ppoints(length(pvector) ))
	plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + stat_abline(intercept=0,slope=1, col="red")
	plot=plot+opts(title=title)
	plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
	plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
	if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
	plot
}

### This is for testing purposes. ######################################
#   set.seed(42)
#   nchr=22
#   nsnps=1000
#   test_data = data.frame(
#       SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
#       CHR=rep(1:nchr,each=nsnps), 
#       BP=rep(1:nsnps,nchr), 
#       P=runif(nchr*nsnps)
#   )
#   ### d[d$SNP=='rs20762',]$P = 1e-29
#   top_snps = c('rs13895','rs20762')
#   surrounding_snps = list(as.character(test_data$SNP[13795:13995]),as.character(test_data$SNP[20662:20862]))
#   
#   pvector = test_data$P
#   names(pvector) = test_data$SNP
#   
#   #CALL:
#   #manhattan(test_data,annotate=top_snps,highlight=surrounding_snps)
#   #qq(pvector,annotate=top_snps,highlight=top_snps)
#########################################################################

# manhattan plot using base graphics
manhattan <- function(dataframe, limitchromosomes=NULL,pt.col=c('gray10','gray50'),pt.bg=c('gray10','gray50'),
    pt.cex=0.45,pch=21,cex.axis=0.7,gridlines=F,gridlines.col='gray83',gridlines.lty=1,gridlines.lwd=1,ymax=8, ymax.soft=T, annotate=NULL,annotate.cex=0.7,annotate.font=3,
    suggestiveline=-log10(1e-5), suggestiveline.col='blue', suggestiveline.lwd=1.5, suggestiveline.lty=1, 
    genomewideline=-log10(5e-8), genomewideline.col='red', genomewideline.lwd=1.5, genomewideline.lty=1, 
    highlight=NULL,highlight.col=c('green3','magenta'),highlight.bg=c('green3','magenta'),  ...) 
{
    #============================================================================================
    ######## Check data and arguments
    d = dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    d = d[( !is.na(d$CHR) & !is.na(d$BP) & !is.na(d$P) ), ]
    
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$CHR)))) warning('non-numeric, non-NA entries in CHR column of dataframe. attempting to remove..')
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$BP)))) warning('non-numeric, non-NA entries in BP column of dataframe. attempting to remove..')
    if (TRUE %in% is.na(suppressWarnings(as.numeric(d$P)))) warning('non-numeric, non-NA entries in P column of dataframe. attempting to remove..')
    
    d = d[!is.na(suppressWarnings(as.numeric(d$CHR))),] # remove rows with non-numeric, non-NA entries
    d = d[!is.na(suppressWarnings(as.numeric(d$BP))),]
    d = d[!is.na(suppressWarnings(as.numeric(d$P))),]
    
    
    if (!is.null(annotate)){
        if ('SNP' %in% names(d)){
                missing_annotate = annotate[!(annotate %in% d$SNP)]
                annotate = annotate[annotate %in% d$SNP]
                
                if (length(missing_annotate)>0){
                        print('These SNPs were not annotated because of missing data:')
                        print(missing_annotate)
                }
        } else {
            stop("D'oh! Dataframe must have a column $SNP with rs_ids to use annotate feature.")
        }
    }
    if (!is.numeric(annotate.cex) | annotate.cex<0) annotate.cex=0.7
    if (!is.numeric(annotate.font)) annotate.font=3
    
    if (is.character(gridlines.col[1]) & !(gridlines.col[1] %in% colors())) gridlines.col = 'gray83'
    if (!is.numeric(pt.cex) | pt.cex<0) pt.cex=0.45
    if (is.character(pt.col) & (FALSE %in% (pt.col %in% colors()))) pt.col = c('gray10','gray50')
    if (is.character(pt.bg) & (FALSE %in% (pt.bg %in% colors()))) pt.bg = F
    if (is.character(highlight.col) & (FALSE %in% (highlight.col %in% colors()))) highlight.col = c('green3','magenta')
    if (is.character(highlight.bg) & (FALSE %in% (highlight.bg %in% colors()))) highlight.bg = F
    if (is.character(suggestiveline.col[1]) & !(suggestiveline.col[1] %in% colors())) suggestiveline.col = 'blue'
    if (is.character(genomewideline.col[1]) & !(genomewideline.col[1] %in% colors())) genomewideline.col = 'red'
        
    if(!is.null(limitchromosomes)){
        if (TRUE %in% is.na(suppressWarnings(as.numeric(limitchromosomes)))){
            stop('limitchromosomes argument is not numeric') 
        } else {  
            d = d[d$CHR %in% as.numeric(limitchromosomes), ]
        }
    }
    

    ######################
    
    # Set positions, ticks, and labels for plotting
    d=subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1)) # sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    
    
    # Ymax
    if(is.na(suppressWarnings(as.numeric(ymax)))){  # not numeric
        ymax = ceiling(max(-log10(d$P)))
        warning('non-numeric ymax argument.')
    } else if (as.numeric(ymax) < 0){           # negative
        ymax = ceiling(max(-log10(d$P)))
        warning('negative ymax argument.')
    }
    if (ymax.soft==T){ #if soft, ymax is just the lower limit for ymax
        ymax = max(ymax, ceiling(max(-log10(d$P))))
        
        # make ymax larger if top annotate SNP is very high
        if (!is.null(annotate)){
            annotate.max = max(d[which(d$SNP %in% annotate),]$logp)
            if ((ymax - annotate.max) < 0.18*ymax){
                ymax = annotate.max + 0.18*ymax
            }
        }
    } #else, ymax = ymax
    
    ## Fix for the bug where one chromosome is missing. Adds index column #####
    d$index=NA
    ind = 0
    for (i in unique(d$CHR)){
        ind = ind + 1
        d[d$CHR==i,]$index = ind
    }
    ########
    
    nchr=length(unique(d$CHR))
    if (nchr==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
        xlabel = paste('Chromosome',unique(d$CHR),'position')
        labs = ticks
    } else {
        ticks = rep(NA,length(unique(d$CHR))+1)
        ticks[1] = 0
        for (i in 1:max(d$index)) {
            d[d$index==i, ]$pos   =    (d[d$index==i, ]$BP - d[d$index==i,]$BP[1]) +1 +ticks[i]
            ticks[i+1] = max(d[d$index==i,]$pos)
        }
        xlabel = 'Chromosome'
        labs = append(unique(d$CHR),'')
    }
    
    # Initialize plot
    xmax = max(d$pos) * 1.03
    xmin = max(d$pos) * -0.03
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03
    plot(0,col=F,xaxt='n',bty='n',xaxs='i',yaxs='i',xlim=c(xmin,xmax), ylim=c(ymin,ymax),
            xlab=xlabel,ylab=expression(-log[10](italic(p))),las=1,cex.axis=cex.axis)
    
    # stagger labels
    blank = rep('',length(labs))
    lowerlabs = rep('',length(labs))
    upperlabs = rep('',length(labs))
    
    for (i in 1:length(labs)){
        if (i %% 2 == 0){
            lowerlabs[i] = labs[i]
        } else{
            upperlabs[i] = labs[i]
        }
    }
    
    axis(1,at=ticks,labels=blank,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
    axis(1,at=ticks,labels=upperlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=-0.25)
    axis(1,at=ticks,labels=lowerlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=0.25)
    
    yvals = par('yaxp')
    yinterval = par('yaxp')[2] / par('yaxp')[3]
    axis(2,at= (seq(0,(ymax+yinterval/2),yinterval) - yinterval/2),labels=F,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
    
    # Gridlines
    if (isTRUE(gridlines)){
        
        abline(v=ticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) #at ticks
        abline(h=seq(0,ymax,yinterval),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) # at labeled ticks
        #abline(h=(seq(0,ymax,yinterval) - yinterval/2),col=gridlines.col[1],lwd=1.0) # at unlabeled ticks
    }
    
    # Points, with optional highlighting
    pt.col = rep(pt.col,max(d$CHR))[1:max(d$CHR)]
    pt.bg = rep(pt.bg,max(d$CHR))[1:max(d$CHR)]
    d.plain = d
    if (!is.null(highlight)) {
        if(class(highlight)!='character' & class(highlight)!='list'){
            stop('"highlight" must be a char vector (for 1 color) or list (for multi color).')
        }
        
        
        if ('SNP' %in% names(d)){
            missing_highlight = highlight[!(highlight %in% d$SNP)]
                highlight = highlight[highlight %in% d$SNP]
                
                if (length(missing_highlight)>0){
                        print('These SNPs were not highlightd because of missing data:')
                        print(missing_highlight)
                }
        
        } else {
            stop("D'oh! Dataframe must have a column $SNP with rs_ids to use highlight feature.")
        }
        
        if (class(highlight)=='character'){ #if char vector, make list for consistency in plotting below
            highlight = list(highlight)
        }
        
        highlight.col = rep(highlight.col,length(highlight))[1:length(highlight)]
        highlight.bg = rep(highlight.bg,length(highlight))[1:length(highlight)]
        
        for (i in 1:length(highlight)){
            d.plain = d.plain[which(!(d.plain$SNP %in% highlight[[i]])), ]
        }
    }
    
    icol=1
    for (i in unique(d.plain$CHR)) {
        with(d.plain[d.plain$CHR==i, ],points(pos, logp, col=pt.col[icol],bg=pt.bg[icol],cex=pt.cex,pch=pch,...))
        icol=icol+1
    }
    
    if (!is.null(highlight)){   
        for (i in 1:length(highlight)){
            d.highlight=d[which(d$SNP %in% highlight[[i]]), ]
            with(d.highlight, points(pos, logp, col=highlight.col[i],bg=highlight.bg[i],cex=pt.cex,pch=pch,...)) 
        }
    }
    
    # Significance lines
    if (is.numeric(suggestiveline)) abline(h=suggestiveline, col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
    if (is.numeric(genomewideline)) abline(h=genomewideline, col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)

    # Annotate
    if (!is.null(annotate)){
        d.annotate = d[which(d$SNP %in% annotate),]
        text(d.annotate$pos,(d.annotate$logp + 0.019*ymax),labels=d.annotate$SNP,srt=90,cex=annotate.cex,adj=c(0,0.48),font=annotate.font)      
    }

    # Box
    box()
}




################################
#### PREPARE DATA FOR GWAS #####
################################



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


# M760.585T667.435 = PC(34:1)
# M758.569T649.183 = PC(34:2)
# M788.616T693.450 = PC(36:1)
# M786.601T676.687== PC(36:2)
# M661.528T579.938 = PE_cer


mPIVUSF <- metabo_p[,c("id","M494.325T326.430","M496.340T365.792","M522.356T384.872","M524.372T425.563","M520.340T347.575","M566.322T352.005","M760.585T667.435","M758.569T649.183","M788.616T693.450","M786.601T676.687","M377.267T388.196","M806.967T347.562","M661.528T579.938","tg","hdl","tc")]

colnames(mPIVUSF)[2:14] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4","PC34_1","PC34_2","PC36_1","PC36_2","MG18_2","LPC18_2bis","PE_cer")


#mPIVUSF$LPC18_1.PC34_1 <-  log2(2^as.numeric(mPIVUSF$LPC18_1)/2^as.numeric(mPIVUSF$PC34_1))
#mPIVUSF$LPC18_1.PC36_1 <-  log2(2^as.numeric(mPIVUSF$LPC18_1)/2^as.numeric(mPIVUSF$PC36_1))


#mPIVUSF$LPC18_0.PC36_1 <-  log2(2^as.numeric(mPIVUSF$LPC18_0)/2^as.numeric(mPIVUSF$PC36_1))
#mPIVUSF$LPC18_0.PC36_2 <-  log2(2^as.numeric(mPIVUSF$LPC18_0)/2^as.numeric(mPIVUSF$PC36_2))

#mPIVUSF$LPC18_2.PC34_2 <-  log2(2^as.numeric(mPIVUSF$LPC18_2)/2^as.numeric(mPIVUSF$PC34_2))
#mPIVUSF$LPC18_2.PC36_2 <-  log2(2^as.numeric(mPIVUSF$LPC18_2)/2^as.numeric(mPIVUSF$PC36_2))


#mPIVUSF$LPC16_0.PC34_1 <-  log2(2^as.numeric(mPIVUSF$LPC18_0)/2^as.numeric(mPIVUSF$PC34_1))
#mPIVUSF$LPC16_0.PC34_2 <-  log2(2^as.numeric(mPIVUSF$LPC18_0)/2^as.numeric(mPIVUSF$PC34_2))


#mPIVUSF$LPC18_2.hdl <-  log2(2^as.numeric(mPIVUSF$LPC18_2)/(mPIVUSF$hdl))
#mPIVUSF$LPC18_2.tc <-  log2(2^as.numeric(mPIVUSF$LPC18_2)/(mPIVUSF$tc))
#mPIVUSF$MG.tg <-  log2(2^as.numeric(mPIVUSF$MG18_2)/mPIVUSF$tg)



# Load .fam file
fam <- read.table("/proj/b2011036/pivus.gwas/pivus.imp.1000gMarch12.b37/pivus-omni-metabo949.b37.sample", stringsAsFactor=F, header=T)
# Cut First line
fam <- fam[2:nrow(fam),1:(ncol(fam)-1)]



pheno <- cbind(paste("PIVUS",sprintf( "%04d", as.numeric(mPIVUSF$id) ),sep=""),
#scale(as.numeric(mPIVUSF$LPC16_0)),
#scale(as.numeric(mPIVUSF$LPC18_0)),
scale(as.numeric(mPIVUSF$LPC18_1)),
scale(as.numeric(mPIVUSF$LPC18_2)),
#scale(as.numeric(mPIVUSF$LPC18_2bis)),
#scale(as.numeric(mPIVUSF$LPC20_4)),
scale(as.numeric(mPIVUSF$MG18_2)),
scale(as.numeric(mPIVUSF$PE_cer)))

#scale(as.numeric(mPIVUSF$LPC18_1.PC34_1)),
#scale(as.numeric(mPIVUSF$LPC18_1.PC36_1)),
#scale(as.numeric(mPIVUSF$LPC18_0.PC36_1)),
#scale(as.numeric(mPIVUSF$LPC18_0.PC36_2)),
#scale(as.numeric(mPIVUSF$LPC18_2.PC34_2)),
#scale(as.numeric(mPIVUSF$LPC18_2.PC36_2)),
#scale(as.numeric(mPIVUSF$LPC16_0.PC34_1)),
#scale(as.numeric(mPIVUSF$LPC16_0.PC34_2)),

#scale(as.numeric(mPIVUSF$LPC18_2.hdl)),
#scale(as.numeric(mPIVUSF$LPC18_2.tc)),
#scale(as.numeric(mPIVUSF$MG.tg)))

#colnames(pheno) <- c("ID_1","LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC18_2bis","LPC20_4","MG18_2","LPC18_1.PC34_1","LPC18_1.PC36_1","LPC18_0.PC36_1","LPC18_0.PC36_2","LPC18_2.PC34_2","LPC18_2.PC36_2","LPC16_0.PC34_1","LPC16_0.PC34_2","LPC18_2.hdl","LPC18_2.tc","MG.tg")

colnames(pheno) <- c("ID_1","LPC18_1","LPC18_2","MG18_2","PE_cer")

phenofam <- merge(fam,pheno, by="ID_1", all.x=T)
phenofam <- apply(phenofam,2,as.character)

phenofam_fin <- phenofam[order(match(phenofam[,1],fam$ID_1)),]

phenofam_finF <- rbind(c(c("0","0","0","D"),rep("P",ncol(pheno)-1)),phenofam_fin)

write.table(phenofam_finF, file="/proj/b2011036/nobackup/pivus_metabo/table4_pivus.sample", quote=F, row.names=F)




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

# M760.585T667.863 = PC(34:1)
# M758.570T649.594 = PC(34:2)
# M788.617T694.035 = PC(36:1)
# M786.601T677.527 = PC(36:2)
# M661.528T580.694 = PE_cer


mTWINGENEF <- metabo_p[,c("twinnr","M494.325T325.856","M496.340T364.064","M522.356T382.722","M524.372T423.532","M520.340T346.671","M566.322T340.659",
"M760.585T667.863","M758.570T649.594","M788.617T694.035","M786.601T677.527",
"M393.231T384.578","M806.968T346.667","M661.528T580.694",colnames(metabo_p)[9756:ncol(metabo_p)])]

colnames(mTWINGENEF)[2:14] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC34_1","PC34_2","PC36_1","PC36_2","MG18_2","LPC18_2bis","PE_cer")



#mTWINGENEF$LPC18_1.PC34_1 <-  log2(2^as.numeric(mTWINGENEF$LPC18_1)/2^as.numeric(mTWINGENEF$PC34_1))
#mTWINGENEF$LPC18_1.PC36_1 <-  log2(2^as.numeric(mTWINGENEF$LPC18_1)/2^as.numeric(mTWINGENEF$PC36_1))


#mTWINGENEF$LPC18_0.PC36_1 <-  log2(2^as.numeric(mTWINGENEF$LPC18_0)/2^as.numeric(mTWINGENEF$PC36_1))
#mTWINGENEF$LPC18_0.PC36_2 <-  log2(2^as.numeric(mTWINGENEF$LPC18_0)/2^as.numeric(mTWINGENEF$PC36_2))

#mTWINGENEF$LPC18_2.PC34_2 <-  log2(2^as.numeric(mTWINGENEF$LPC18_2)/2^as.numeric(mTWINGENEF$PC34_2))
#mTWINGENEF$LPC18_2.PC36_2 <-  log2(2^as.numeric(mTWINGENEF$LPC18_2)/2^as.numeric(mTWINGENEF$PC36_2))


#mTWINGENEF$LPC16_0.PC34_1 <-  log2(2^as.numeric(mTWINGENEF$LPC18_0)/2^as.numeric(mTWINGENEF$PC34_1))
#mTWINGENEF$LPC16_0.PC34_2 <-  log2(2^as.numeric(mTWINGENEF$LPC18_0)/2^as.numeric(mTWINGENEF$PC34_2))


#mTWINGENEF$LPC18_2.hdl <-  log2(2^as.numeric(mTWINGENEF$LPC18_2)/(mTWINGENEF$hdl))
#mTWINGENEF$LPC18_1.hdl <-  log2(2^as.numeric(mTWINGENEF$LPC18_1)/(mTWINGENEF$hdl))
#mTWINGENEF$LPC18_2.tc <-  log2(2^as.numeric(mTWINGENEF$LPC18_2)/(mTWINGENEF$tc))
#mTWINGENEF$MG.tg <-  log2(2^as.numeric(mTWINGENEF$MG18_2)/mTWINGENEF$tg)



fam <- read.table("/proj/b2011036/twge.gwas/twge.imp.1000gMarch12.b37/meta/fam/twge.+.qced.fam", stringsAsFactor=F, header=F)
pca <- read.table("/proj/b2011036/twge.gwas/twge.imp.hapmap.b36/twge.+.qced.mds", stringsAsFactor=F, header=T)

fam2 <- merge(fam,pca, by.x="V2", by.y="IID", all.x=T)


pheno <- cbind(mTWINGENEF$gwas_id,mTWINGENEF$age,
#scale(as.numeric(mTWINGENEF$LPC16_0)),
#scale(as.numeric(mTWINGENEF$LPC18_0)),
scale(as.numeric(mTWINGENEF$LPC18_1)),
scale(as.numeric(mTWINGENEF$LPC18_2)),
#scale(as.numeric(mTWINGENEF$LPC18_2bis)),
#scale(as.numeric(mTWINGENEF$LPC20_4)),
scale(as.numeric(mTWINGENEF$MG18_2)),
scale(as.numeric(mTWINGENEF$PE_cer)))


#scale(as.numeric(mTWINGENEF$LPC18_1.PC34_1)),
#scale(as.numeric(mTWINGENEF$LPC18_1.PC36_1)),
#scale(as.numeric(mTWINGENEF$LPC18_0.PC36_1)),
#scale(as.numeric(mTWINGENEF$LPC18_0.PC36_2)),
#scale(as.numeric(mTWINGENEF$LPC18_2.PC34_2)),
#scale(as.numeric(mTWINGENEF$LPC18_2.PC36_2)),
#scale(as.numeric(mTWINGENEF$LPC16_0.PC34_1)),
#scale(as.numeric(mTWINGENEF$LPC16_0.PC34_2)),

#scale(as.numeric(mTWINGENEF$LPC18_2.hdl)),
#scale(as.numeric(mTWINGENEF$LPC18_1.hdl)),
#scale(as.numeric(mTWINGENEF$LPC18_2.tc)),
#scale(as.numeric(mTWINGENEF$MG.tg)))

#colnames(pheno) <- c("V2","Age","LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC18_2bis","LPC20_4","MG18_2","LPC18_1.PC34_1","LPC18_1.PC36_1","LPC18_0.PC36_1","LPC18_0.PC36_2","LPC18_2.PC34_2","LPC18_2.PC36_2","LPC16_0.PC34_1","LPC16_0.PC34_2","LPC18_2.hdl","LPC18_2.tc","MG.tg")

colnames(pheno) <- c("V2","Age","LPC18_1","LPC18_2","MG18_2","PE_cer")


phenofam <- merge(fam2,pheno, by="V2", all.x=T)
phenofam <- apply(phenofam,2,as.character)

# Missin = -99
phenofam[is.na(phenofam)] <- -99

# Keep right columns and reorder
phenofam_fin <- phenofam[,c(2,1,5,9,10,11,13:ncol(phenofam))]
colnames(phenofam_fin)[1:7] <- c("FID","IID","Gender","C1","C2","C3","Age")

# Order as the original fam file
phenofam_finF <- phenofam_fin[order(match(phenofam_fin[,2],fam$V2)),]



# Write sample file for plink
write.table(phenofam_finF, file="/proj/b2011036/nobackup/twge_metabo/table4_twge.sample", quote=F, row.names=F)



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


# M760.594T674.770 = PC(34:1)
# M758.569T657.575 = PC(34:2)
# M788.616T700.812 = PC(36:1)
# M786.609T683.883 = PC(36:2)
# M661.528T590.990 = PE_cer

mULSAMF <- metabo_p[metabo_p$time==0 & metabo_p$double_==0,c("pat_time","M494.324T324.997","M496.340T375.465","M522.356T395.236","M524.372T436.589","M520.340T357.723","M566.322T351.269",
"M760.594T674.770","M758.569T657.575","M788.616T700.812","M786.609T683.883",
"M393.238T397.936","M806.966T357.935","M661.528T590.990",colnames(metabo_p)[10163:ncol(metabo_p)])]

colnames(mULSAMF)[2:14] <- c("LPC16_1","LPC16_0","LPC18_1","LPC18_0","LPC18_2","LPC20_4",
"PC34_1","PC34_2","PC36_1","PC36_2","MG18_2","LPC18_2bis","PE_cer")


## Make ratios 

#mULSAMF$LPC16_0.PC34_1 <-  log2(2^as.numeric(mULSAMF$LPC18_0)/2^as.numeric(mULSAMF$PC34_1))
#mULSAMF$LPC16_0.PC34_2 <-  log2(2^as.numeric(mULSAMF$LPC18_0)/2^as.numeric(mULSAMF$PC34_2))

#mULSAMF$LPC18_1.PC34_1 <-  log2(2^as.numeric(mULSAMF$LPC18_1)/2^as.numeric(mULSAMF$PC34_1))
#mULSAMF$LPC18_1.PC36_1 <-  log2(2^as.numeric(mULSAMF$LPC18_1)/2^as.numeric(mULSAMF$PC36_1))

#mULSAMF$LPC18_2.PC34_2 <-  log2(2^as.numeric(mULSAMF$LPC18_2)/2^as.numeric(mULSAMF$PC34_2))
#mULSAMF$LPC18_2.PC36_2 <-  log2(2^as.numeric(mULSAMF$LPC18_2)/2^as.numeric(mULSAMF$PC36_2))

#mULSAMF$LPC18_0.PC36_1 <-  log2(2^as.numeric(mULSAMF$LPC18_0)/2^as.numeric(mULSAMF$PC36_1))
#mULSAMF$LPC18_0.PC36_2 <-  log2(2^as.numeric(mULSAMF$LPC18_0)/2^as.numeric(mULSAMF$PC36_2))

#mULSAMF$LPC18_2.hdl <-  log2(2^as.numeric(mULSAMF$LPC18_2)/(mULSAMF$hdl))
#mULSAMF$LPC18_1.hdl <-  log2(2^as.numeric(mULSAMF$LPC18_1)/(mULSAMF$hdl))
#mULSAMF$LPC18_2.tc <-  log2(2^as.numeric(mULSAMF$LPC18_2)/(mULSAMF$tc))
#mULSAMF$MG.tg <-  log2(2^as.numeric(mULSAMF$MG18_2)/mULSAMF$tg)




# Load .fam file
fam <- read.table("/proj/b2011036/ulsam.gwas/ulsam.imp.1000gMarch12.b37/ULSAMb37.sample", stringsAsFactor=F, header=T)
# Cut First line
fam <- fam[2:nrow(fam),1:(ncol(fam)-1)]


pheno <- cbind(paste("ULSAM",sprintf( "%04d", as.numeric(mULSAMF$pat) ),sep=""),
#scale(as.numeric(mULSAMF$LPC16_0)),
#scale(as.numeric(mULSAMF$LPC18_0)),
scale(as.numeric(mULSAMF$LPC18_1)),
scale(as.numeric(mULSAMF$LPC18_2)),
#scale(as.numeric(mULSAMF$LPC18_2bis)),
#scale(as.numeric(mULSAMF$LPC20_4)),
scale(as.numeric(mULSAMF$MG18_2)),
scale(as.numeric(mULSAMF$PE_cer)))

#scale(as.numeric(mULSAMF$LPC18_1.PC34_1)),
#scale(as.numeric(mULSAMF$LPC18_1.PC36_1)),
#scale(as.numeric(mULSAMF$LPC18_0.PC36_1)),
#scale(as.numeric(mULSAMF$LPC18_0.PC36_2)),
#scale(as.numeric(mULSAMF$LPC18_2.PC34_2)),
#scale(as.numeric(mULSAMF$LPC18_2.PC36_2)),
#scale(as.numeric(mULSAMF$LPC16_0.PC34_1)),
#scale(as.numeric(mULSAMF$LPC16_0.PC34_2)),

#scale(as.numeric(mULSAMF$LPC18_2.hdl)),
#scale(as.numeric(mULSAMF$LPC18_1.hdl)))
#scale(as.numeric(mULSAMF$LPC18_2.tc)),
#scale(as.numeric(mULSAMF$MG.tg)))

#colnames(pheno) <- c("ID_1","LPC16_0","LPC18_0","LPC18_1","LPC18_2","LPC18_2bis","LPC20_4","MG18_2","LPC18_1.PC34_1","LPC18_1.PC36_1","LPC18_0.PC36_1","LPC18_0.PC36_2","LPC18_2.PC34_2","LPC18_2.PC36_2","LPC16_0.PC34_1","LPC16_0.PC34_2","LPC18_2.hdl","LPC18_2.tc","MG.tg")

colnames(pheno) <- c("ID_1","LPC18_1","LPC18_2","MG18_2","PE_cer")

phenofam <- merge(fam,pheno, by="ID_1", all.x=T)
phenofam <- apply(phenofam,2,as.character)

phenofam_fin <- phenofam[order(match(phenofam[,1],fam$ID_1)),]

phenofam_finF <- rbind(c(c("0","0","0","D"),rep("P",ncol(pheno)-1)),phenofam_fin)


write.table(phenofam_finF, file="/proj/b2011036/nobackup/ulsam_metabo/table4_ulsam.sample", quote=F, row.names=F)



##########################
##### READ GWAS DATA #####
##########################



### PHENOTYPE AND SPECIFIC SNPS ###
phs=c("LPC18_1","LPC18_2","MG18_2","PE_cer")



#### PIVUS ####
for (ph in phs)
{
	cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), header=T, nrows=5)
	classes <- sapply(cc, class)
	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), intern=T))," ")[[1]][1]
	t <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes)

	t2 <- t[t[,12] > 0 & t[,12] < 1 & t[,14] > 0 & t[,6] > 0.4 & t[,10] >= 0.03 & t[,11]>1e-6, ] 
	t2 <- cbind(paste("chr",t2[,1],":",t2[,3],sep=""),t2)
	
	t2$freq <- (t2$all_BB+1/2*t2$all_AB)/(t2$all_AA+t2$all_AB+t2$all_BB)
	
	write.csv(t2[t2[,13]<5e-04,], file=paste("/proj/b2011036/nobackup/",ph,".final.pivus_e_4.csv", sep=""))

	# Save for metal
	t3 <- t2[,c(1,5,6,16,13,14,15)]
	colnames(t3) <- c("markername","other_allele","effect_allele","eaf","p","beta","se")
	write.table(t3, file=paste("/proj/b2011036/nobackup/",ph,".final.pivus_metal", sep=""),col.names=T, row.names=F, quote=F)
	
	print(ph)
}





#### ULSAM ####

for (ph in phs)
{
	cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), header=T, nrows=5)
	classes <- sapply(cc, class)
	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), intern=T))," ")[[1]][1]
	t <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows))

	t2 <- t[t[,12] > 0 & t[,12] < 1 & t[,14] > 0 & t[,6] > 0.4 & t[,10] >= 0.03 & t[,11]>1e-6, ] 
	t2 <- cbind(paste("chr",t2[,1],":",t2[,3],sep=""),t2)
	
	t2$freq <- (t2$all_BB+1/2*t2$all_AB)/(t2$all_AA+t2$all_AB+t2$all_BB)
	
	
	write.csv(t2[t2[,13]<5e-04,], file=paste("/proj/b2011036/nobackup/",ph,".final.ulsam_e_4.csv", sep=""))

	# Save for metal
	t3 <- t2[,c(1,5,6,16,13,14,15)]
	colnames(t3) <- c("markername","other_allele","effect_allele","eaf","p","beta","se")
	write.table(t3, file=paste("/proj/b2011036/nobackup/",ph,".final.ulsam_metal", sep=""),col.names=T, row.names=F, quote=F)
	
	
	print(ph)
}





#### TWINGENE ####



for (ph in phs)
{
	cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), header=T, nrows=5)
	classes <- sapply(cc, class)
	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), intern=T))," ")[[1]][1]
	t <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes)

	t2 <- t[!is.na(t[,10]) & t[,10]>0 & t[,10]<1 & t[,9]>0 & t[,7] > 0.4 & (t[,6] >= 0.03 & t[,6] <= 0.97), ] 
	t2 <- cbind(paste("chr",t2[,1],":",t2[,3],sep=""),t2)
	
	write.csv(t2[t2[,11]<5e-04,], file=paste("/proj/b2011036/nobackup/",ph,".final.twge_e_4.csv", sep=""))

	# Save for metal
	t3 <- t2[,c(1,5,6,7,11,9,10)]
	colnames(t3) <- c("markername","effect_allele","other_allele","eaf","p","beta","se")
	write.table(t3, file=paste("/proj/b2011036/nobackup/",ph,".final.twge_metal", sep=""),col.names=T, row.names=F, quote=F)
	
	print(ph)
}

#####################################
##### READ META-ANALYSIS RESULTS ####
#####################################


phs=c("LPC18_1","LPC18_2","MG18_2","PE_cer")



for (ph in phs)
{
	cc <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, nrows=5)
	classes <- sapply(cc, class)
	nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), intern=T))," ")[[1]][1]
	t <- read.table(paste("/proj/b2011036/nobackup/",ph,"_metaanalysis1.tbl", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes, stringsAsFactor=F)

	write.csv(t[t[,10]<5e-06,], file=paste("/proj/b2011036/nobackup/",ph,".metaanalysis_e_6.csv", sep=""))
	
	#spl <- strsplit(as.character(t$MarkerName),":")
	
	#forman <- data.frame(SNP=as.character(t$MarkerName),CHR=as.numeric(gsub("chr","",sapply(spl,"[[",1))), BP=as.numeric(sapply(spl,"[[",2)), P=t$P.value)	
	#bitmap(paste("/proj/b2011036/nobackup/",ph,".metaanalysis_manhattan.bmp", sep=""), res=300)
	#manhattan(forman,pt.col=c("black","#666666","#CC6600"), pch=20)
	#dev.off()
	#print(ph)
	
}







### Top-SNP for LysoPC 18:1 has significant P-value for heterogenity, try random-effect meta-analysis ###

ph <- "MG18_2"

cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), intern=T))," ")[[1]][1]
tP <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.pivus.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes)




cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), intern=T))," ")[[1]][1]
tU <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.ulsam.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows))

cc <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), header=T, nrows=5)
classes <- sapply(cc, class)
nrows=strsplit(as.character(system(paste("wc -l /proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), intern=T))," ")[[1]][1]
tT <- read.table(paste("/proj/b2011036/nobackup/",ph,".final.twge.txt", sep=""), header=T, comment.char = "", nrow=as.numeric(nrows), colClasses=classes)


## chr8:94088655

beta <- c(tP[tP[,1]==8 & tP[,3]==94088655,13],tU[tU[,1]==8 & tU[,3]==94088655,13],-tT[tT[,1]==8 & tT[,3]==94088655,8])
se <- c(tP[tP[,1]==8 & tP[,3]==94088655,14],tU[tU[,1]==8 & tU[,3]==94088655,14],tT[tT[,1]==8 & tT[,3]==94088655,9])

library(meta)
metaana <- metagen(beta,se,studlab=c("Pivus","Ulsam","TwinGene"))








