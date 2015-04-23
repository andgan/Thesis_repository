#----------------------------------------------
# Filename: code_variables_cat.R
# Study: UK mortality
# Author: Andrea Ganna
# Date: 13JUN2014
# Updated: 
# Purpose: clean the initial data and create a version of the dataset to be used for imputation. It also tranform all the continuous variable in quantiles with constrain.
# Note: 
#-----------------------------------------------
# Data used: out4.Rdata, var_and_codes.csv 
# Data created: out5.Rdata and Descriptives_statistics.pdf
#-----------------------------------------------
# OP: R 3.1.0
#-----------------------------------------------*/


## Load in-house functions
source("/proj/b2011036/uk.biobank/Pgm/new_function.R")

## Load functions
load("/proj/b2011036/uk.biobank/pre_imputation/out4.Rdata")

# Load variables names (used in the final plot) #
ph <- read.csv("/proj/b2011036/uk.biobank/var_and_codes.csv", header=T, stringsAsFactor=F)

### EXTERNAL VARIABLES ####
bdE4$center <- bdE4$f.54.0.0
bdE4$date_center <- bdE4$f.53.0.0
bdE4$sex <- bdE4$f.31.0.0
bdE4$age <- bdE4$f.21003.0.0
bdE4$ddate <- bdE4$f.40000.0.0
bdE4$out <- ifelse(is.na(bdE4$ddate),0,1)
# Cause of deaths, code in 5 categories
t1 <- as.numeric(substr(bdE4$f.40001.0.0,2,3))
t2 <- substr(bdE4$f.40001.0.0,0,1)
bdE4$cofd<- as.factor(ifelse(t2=="C" | ( t2=="D" & t1 < 49 ),"Cancer",
ifelse(t2=="I" ,"Cardiovascular",
ifelse(t2=="J", "Respiratory",
ifelse(t2=="K", "Digestive",
ifelse(t2=="V" | t2=="W" | t2=="X" | t2=="Y","External causes",
"other diseases"))))))


## CHARLSON's SCORE - SELF-REPORTED ##
c1 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))==1075}),1,any) # heart attack/myocardial infarction
c1[is.na(c1)] <- 0

c2 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1076,1079,1588)}),1,any) # heart failure/pulmonary odema; cardiomyopathy: hypertrophic cardiomyopathy (hcm / hocm)
c2[is.na(c2)] <- 0

c3 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1067,1087,1088, 1492,1591,1592)}),1,any) # peripheral vascular disease: leg claudication/ intermittent claudication, arterial embolism; aortic aneurysm: aortic aneurysm rupture, aortic dissection
c3[is.na(c3)] <- 0

c4 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1081,1583,1491,1086,1082,1083,1425)}),1,any) # cerebrovascular disease: stroke: ischaemic stroke, brain haemorrhage, subarachnoid haemorrhage; transient ischaemic attack (tia), subdural haemorrhage/haematoma, cerebral aneurysm
c4[is.na(c4)] <- 0

c5 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))==1263}),1,any) # dementia/alzheimers/cognitive impairment
c5[is.na(c5)] <- 0

c6 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1112,1111,1113,1412,1472,1114)}),1,any) # chronic obstructive airways disease/copd, asthma; emphysema/chronic bronchitis: bronchitis, emphysema; bronchiectasis
c6[is.na(c6)] <- 0

c7 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1464, 1381,1383,1480,1481,1384, 1372,1376, 1377, 1378, 1379, 1380, 1373)}),1,any) #rheumatoid arthritis, systemic lupus erythematosis/sle, dermatopolymyositis: dermatomyositis, polymyositis; scleroderma/systemic sclerosis; vasculitis: giant cell/temporal arteritis, polymyalgia rheumatica, wegners granulmatosis, microscopic polyarteritis,polyartertis nodosa; connective tissue disorder
c7[is.na(c7)] <- 0

c8 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1400,1142,1457)}),1,any) # peptic ulcer, gastric/stomach ulcers, duodenal ulcer
c8[is.na(c8)] <- 0

c9 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1155, 1156, 1157, 1578, 1579, 1580, 1581, 1582, 1136, 1158, 1506,1604 )}),1,any) # liver disease: hepatitis: infective/viral hepatitis, hepatitis a, hepatitis b, hepatitis c, hepatitis d, hepatitis e, non-infective hepatitis. liver failure/cirrhosis: primary biliary cirrhosis, alcoholic liver disease / alcoholic cirrhosis; liver/biliary/pancreas problem.
c9[is.na(c9)] <- 0

c10 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1220,1222,1223)}),1,any) #Diabetes: type 1 diabetes, type 2 diabetes
c10[is.na(c10)] <- 0

c11 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1276,1468,1607)}),1,any) #Diabetes with end-organ damage: diabetic eye disease, diabetic neuropathy/ulcers, diabetic nephropathy 
c11[is.na(c11)] <- 0
c11[c11==1] <- 2

c12 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))==c(1252)}),1,any) #Paraplegia
c12[is.na(c12)] <- 0
c12[c12==1] <- 2

c13 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))%in%c(1192,1193,1194,1608,1609)}),1,any) # renal/kidney failure: renal failure requiring dialysis, renal failure not requiring dialysis; nephritis: glomerulnephritis
c13[is.na(c13)] <- 0
c13[c13==1] <- 2

c14 <- apply(sapply(0:5,function(i) {!eval(parse(text=paste0("bdE4$f.20001.0.",i)))%in%c(NA,1060 ,1061,1073,1062,1071,1085)}),1,any) # solid tumor and exclude non-melanoma skin cancer: basal cell carcinoma, rodent ulcer, squamous cell carcinoma; metastatic cancer (unknown primary): bone metastases / bony secondaries; 
c14[is.na(c14)] <- 0
c14[c14==1] <- 2

c15 <- apply(sapply(0:5,function(i) {eval(parse(text=paste0("bdE4$f.20001.0.",i)))%in%c(1071,1085)}),1,any) #metastatic cancer (unknown primary): bone metastases / bony secondaries; 
c15[is.na(c15)] <- 0
c15[c15==1] <- 3

c16 <- apply(sapply(0:28,function(i) {eval(parse(text=paste0("bdE4$f.20002.0.",i)))==1439}),1,any) #hiv/aids
c16[is.na(c16)] <- 0
c16[c16==1] <- 6


bdE4$charlsonSR <- c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13+c14+c15+c16
#write.csv(cbind(bdE4$f.eid,charlsonSR), file="Charlson.csv")

## CHARLSON's SCORE - REGISTRIES ##

#h1 <- apply(sapply(0:263, function(i) {grepl(paste(c('I21','I22','I252'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Myocardial Infarction

#h2 <- apply(sapply(0:263, function(i) {grepl(paste(c('I43','I50','I099','I110','I130','I132',
#'I255','I420','I425','I426','I427','I428','I429','P290'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Congestive Heart Failure

#h3 <- apply(sapply(0:263, function(i) {grepl(paste(c('I70','I71', 'I731','I738','I739','I771',
#'I790','I792','K551','K558','K559','Z958','Z959'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Periphral Vascular Disease

#h4 <- apply(sapply(0:263, function(i) {grepl(paste(c('G45','G46','I60','I61','I62','I63','I64','I65','I66','I67','I68','I69','H340'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Cerebrovascular Disease

#h5 <- apply(sapply(0:263, function(i) {grepl(paste(c('F00','F01','F02','F03','G30','F051','G311'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Dementia

#h6 <- apply(sapply(0:263, function(i) {grepl(paste(c('J40','J41','J42','J43','J44','J45','J46','J47',
#'J60','J61','J62','J63','J64','J65','J66','J67',
#'I278','I279','J684','J701','J703'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Chronic Pulmonary Disease

#h7 <- apply(sapply(0:263, function(i) {grepl(paste(c('M05','M32','M33','M34','M06','M315','M351','M353','M360'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Connective Tissue Disease-Rheumatic Disease

#h8 <- apply(sapply(0:263, function(i) {grepl(paste(c('K25','K26','K27','K28'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Peptic Ulcer Disease 

#h9 <- apply(sapply(0:263, function(i) {grepl(paste(c('B18','K73','K74','K700','K701','K702','K703','K709',
#'K717','K713','K714','K715','K760','K762','K763','K764','K768','K769','Z944'),collapse="|"), #eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Mild Liver Disease

#h10 <- apply(sapply(0:263, function(i) {grepl(paste(c('E100','E101','E106','E108','E109','E110','E111','E116','E118','E119',
#'E120','E121','E126','E128','E129',
#'E130','E131','E136','E138','E139',
#'E140','E141','E146','E148','E149'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any)+0 # Diabetes without complications

#h11 <- apply(sapply(0:263, function(i) {grepl(paste(c('E102','E103','E104','E105','E107',
#'E112','E113','E114','E115','E117',
#'E122','E123','E124','E125','E127',
#'E132','E133','E134','E135','E137',
#'E142','E143','E144','E145','E147'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Diabetes with complications
#h11 <- ifelse(h11==T,2,0)

#h12 <- apply(sapply(0:263, function(i) {grepl(paste(c('G81','G82','G041','G114','G801','G802',
#'G830','G831','G832','G833','G834','G839'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Paraplegia and Hemiplegia
#h12 <- ifelse(h12==T,2,0)
 
#h13 <- apply(sapply(0:263, function(i) {grepl(paste(c('N18','N19','N052','N053','N054','N055','N056','N057',
#'N250','I120','I131','N032','N033','N034','N035','N036','N037',
#'Z490','Z491','Z492','Z940','Z992'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Renal Disease
#h13 <- ifelse(h13==T,2,0)

#h14 <- apply(sapply(0:263, function(i) {grepl(paste(c('C00','C01','C02','C03','C04','C05','C06','C07','C08','C09',
#'C10','C11','C12','C13','C14','C15','C16','C17','C18','C19',
#'C20','C21','C22','C23','C24','C25','C26',
#'C30','C31','C32','C33','C34','C37','C38','C39',
#'C40','C41','C43','C45','C46','C47','C48','C49',
#'C50','C51','C52','C53','C54','C55','C56','C57','C58',
#'C60','C61','C62','C63','C64','C65','C66','C67','C68','C69',
#'C70','C71','C72','C73','C74','C75','C76',
#'C81','C82','C83','C84','C85','C88',
#'C90','C91','C92','C93','C94','C95','C96','C97'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Cancer
#h14 <- ifelse(h14==T,2,0)

#h15 <- apply(sapply(0:263, function(i) {grepl(paste(c('K704','K711','K721','K729','K765','K766','K767','I850','I859','I864','I982'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Moderate or Severe Liver Disease
#h15 <- ifelse(h15==T,3,0)

#h16 <- apply(sapply(0:263, function(i) {grepl(paste(c('C77','C78','C79','C80'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # Metastatic Carcinoma
#h16 <- ifelse(h16==T,3,0)

#h17 <- apply(sapply(0:263, function(i) {grepl(paste(c('B20','B21','B22','B24'),collapse="|"), eval(parse(text=paste0("bdE4$f.41202.0.",i))))}),1,any) # AIDS/HIV
#h17 <- ifelse(h17==T,6,0)

#charlsonREG <- h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12+h13+h14+h15+h16+h17

#summary(glm(bdE4$out~bdE4$age+bdE4$sex+charlsonREG, family=binomial()))
#summary(glm(bdE4$out~bdE4$age+bdE4$sex+charlsonSR, family=binomial()))


## Redefine cut function to include only unique cuts and include lowest cut ##
cut2 <- function(var,out=bdE4$out,min_op_r=0.20,min_tg_n=20){	
	optB <- optimal_bin(var,out,min_op_r,min_tg_n)
	return(cut(var,unique(c(min(optB[,2]),optB[,3])),include.lowest = T))
	}


##########################################
##### CHECK EACH VARIABLE AND PROCESS ####
##########################################

# THE FOLLOWING VARIABLES SHOULD BE CODED LATER (LARGER CODE), BUT I DO IT HERE BECAUSE THEY ARE USEFUL LATER #

# Qualifications #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6138" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6138.0.0)))
	{
		ll <- levels(bdE4$f.6138.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6138.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6138.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6138" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Gas or solid−fuel cooking/heating #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6139" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6139.0.0)))
	{
		ll <- levels(bdE4$f.6139.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6139.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6139.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6139" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Heating type(s) in home #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6140" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6140.0.0)))
	{
		ll <- levels(bdE4$f.6140.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6140.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6140.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6140" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# How are people in household related to participant #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6141" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No house or live alone")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.670.0.0 %in% c("Sheltered accommodation","Care home") | bdE4$f.709.0.0==1] <- "No house or live alone"
		bdE4[paste("f.6141.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6141.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6141" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Current employment status #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6142" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6142.0.0)))
	{
		ll <- levels(bdE4$f.6142.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6142.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6142.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6142" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Transport type for commuting to job workplace #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6143" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{	
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No working")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.6142.0.l2 == 0] <- "No working"
		bdE4[paste("f.6143.0.l",k,sep="")] <- tef	
		levels(bdE4[[paste("f.6143.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6143" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Never eat eggs, dairy, wheat, sugar #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6144" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6144.0.0)))
	{
		ll <- levels(bdE4$f.6144.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6144.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6144.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6144" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Illness, injury, bereavement, stress in last 2 years #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6145" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6145.0.0)))
	{
		ll <- levels(bdE4$f.6145.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6145.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6145.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6145" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Attendance/disability/mobility allowance #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6146" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6146.0.0)))
	{
		ll <- levels(bdE4$f.6146.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6146.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6146.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6146" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Reason for glasses/contact lenses #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6147" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No wearing glasses")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.2207.0.0 != "Yes" ] <- "No wearing glasses"
		bdE4[paste("f.6147.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6147.0.l",k,sep="")]])[2] <- ll	
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6147" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Eye problems/disorders #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6148" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6148.0.0)))
	{
		ll <- levels(bdE4$f.6148.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6148.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6148.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6148" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Mouth/teeth dental problems #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6149" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6149.0.0)))
	{
		ll <- levels(bdE4$f.6149.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6149.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6149.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6149" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Vascular/heart problems diagnosed by doctor #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6150" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6150.0.0)))
	{
		ll <- levels(bdE4$f.6150.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6150.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6150.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6150" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Fractured bone site(s) #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6151" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6151.0.0)))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No fractured bone")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.2463.0.0 != "Yes"] <- "No fractured bone"
		bdE4[paste("f.6151.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6151.0.l",k,sep="")]])[2] <- ll	
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6151" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# DVT, bronchitis, emphysema, asthma etc.. #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6152" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6152.0.0)))
	{
		ll <- levels(bdE4$f.6152.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6152.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6152.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6152" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Blood pressure meds - cholesterol etc. #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6153" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6153.0.0)))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$sex != "Female" ] <- NA
		bdE4[paste("f.6153.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6153.0.l",k,sep="")]])[2] <- ll	
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6153" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Medication for pain relief, constipation, heartburn #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6154" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6154.0.0)))
	{
		ll <- levels(bdE4$f.6154.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6154.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6154.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6154" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Vitamin and mineral suppl #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6155" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6155.0.0)))
	{
		ll <- levels(bdE4$f.6155.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6155.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6155.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6155" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Manic/hyper symptoms #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6156" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"I'm not manic/irritable")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.4642.0.0 != "Yes"  & bdE4$f.4653.0.0!="Yes"] <- "I'm not manic/irritable"
		bdE4[paste("f.6156.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6156.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6156" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Why stopped smoking #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6157" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No past smokers or stop too short")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.1249.0.0 != "Smoked on most or all days"  | bdE4$f.2907.0.0!="Yes"] <- "No past smokers or stop too short"
		bdE4[paste("f.6157.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6157.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6157" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]


# Why reduced smoking #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6158" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No smokers or stop too short")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.1239.0.0 != "Yes, on most or all days"  | bdE4$f.3506.0.0!="Less nowadays?"] <- "No smokers or stop too short"
		bdE4[paste("f.6158.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6158.0.l",k,sep="")]])[2] <- ll	
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6158" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]


# Pain type(s) experienced in last month #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6159" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6159.0.0)))
	{
		ll <- levels(bdE4$f.6159.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6159.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6159.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6159" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Leisure/social activities #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6160" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6160.0.0)))
	{
		ll <- levels(bdE4$f.6160.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6160.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6160.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6160" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Types of transport used #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6162" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"No walking")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.864.0.0 == -2 ] <- "No walking"
		bdE4[paste("f.6162.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6162.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6162" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]


# Types of physical activity in last 4 weeks #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6164" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6164.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6164.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6164" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Medication for cholesterol, blood pressure or diabetes #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6177" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6177.0.0)))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$sex != "Male" ] <- NA
		bdE4[paste("f.6177.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.6177.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6177" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

# Mineral and other dietary supplements #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="6179" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(bdE4$f.6179.0.0)))
	{
		ll <- levels(bdE4$f.6179.0.0)[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		te[rowSums(is.na(tt))==ncol(tt)] <- NA
		bdE4[paste("f.6179.0.l",k,sep="")] <- as.factor(te*1)
		levels(bdE4[[paste("f.6179.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="6179" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]

print("first")

### Time taken to different tests ### 
bdE4$f.3.0.0 <- NULL
bdE4$f.4.0.0 <- NULL
bdE4$f.5.0.0 <- NULL
bdE4$f.6.0.0 <- NULL

## Heel ultrasound method #
bdE4$f.19.0.0 <- NULL

# Weight method #
bdE4$f.21.0.0 <- NULL

# Spirometry method #
bdE4$f.23.0.0 <- NULL

# Sex - coded as external variable #
bdE4$f.31.0.0 <- NULL

# year of birth #
bdE4$f.34.0.0 <- NULL

# Blood pressure device ID #
bdE4$f.36.0.0 <- NULL

# Blood pressure manual sphygmomanometer device ID #
bdE4$f.37.0.0 <- NULL

# Month of birth #
bdE4$f.52.0.0 <- NULL

# Date attending assessment center - coded as external variable #
bdE4$f.53.0.0  <- NULL

# Assessment center - coded as external variable #
bdE4$f.54.0.0  <- NULL

# Month Assessment center #
bdE4$f.55.0.0  <- NULL

# Willing to attempt cognitive tests #
bdE4$f.62.0.0  <- NULL

# EXCLUSION SEVERAL VARIABLES WITH YEAR/AGE FOR CANCER #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("84","87","92")]

# Pulse rate during blood pressure measurement, not useful and small sample size #
bdE4$f.95.0.0 <- NULL
bdE4$f.95.0.1 <- NULL

# Time form interview start to blood pressure is measured #
bdE4$f.96.0.0 <- NULL
bdE4$f.96.0.1 <- NULL

# Pulse rate - automatic reading - average between two measurements #
bdE4$f.102.0.0 <- rowMeans(cbind(bdE4$f.102.0.0,bdE4$f.102.0.1), na.rm=T)
bdE4$f.102.0.1 <- NULL

# Aide−memoire completed #
bdE4$f.111.0.0 <- NULL

# Birth weight known #
bdE4$f.120.0.0 <- as.factor(ifelse(bdE4$f.120.0.0=="Yes - pounds and ounces" | bdE4$f.120.0.0=="Yes - Kilograms","Yes","No"))

# Code Jobs in 9 main categories - use job code deduced when free text entry # 
t <- substr(bdE4$f.132.0.0,0,1)
t <- ifelse(t==0,substr(bdE4$f.20024.0.0,0,1),t)
bdE4$f.132.0.M <- as.factor(ifelse(t==1,"Managers and Senior Officials",
ifelse(t==2,"Professional Occupations",
ifelse(t==3,"Associate Professional and Technical Occupations",
ifelse(t==4,"Administrative and Secretarial Occupations",
ifelse(t==5,"Skilled Trades Occupations",
ifelse(t==6,"Personal Service Occupations",
ifelse(t==7,"Sales and Customer Service Occupations",
ifelse(t==8,"Process, Plant and Machine Operatives",
ifelse(t==9,"Elementary Occupations",t))))))))))
levels(bdE4$f.132.0.M) <- c(levels(bdE4$f.132.0.M),"Not working")
bdE4$f.132.0.M[bdE4$f.6142.0.l2==0] <- "Not working"
bdE4$f.132.0.0 <- NULL
bdE4$f.20024.0.0 <- NULL

# N. of self reported cancers or illnesses or operations #
bdE4$f.134.0.0b <- bdE4$f.134.0.0
bdE4$f.135.0.0b <- bdE4$f.135.0.0
bdE4$f.136.0.0b <- bdE4$f.136.0.0
bdE4$f.137.0.0b <- bdE4$f.137.0.0

# Make this in categories and factorial otherwise the cutting does not work #
bdE4$f.134.0.0 <- as.factor(ifelse(bdE4$f.134.0.0>=2,">2",bdE4$f.134.0.0))

# N. of culumn shown in the pair game #
bdE4$f.396.0.1 <- NULL
bdE4$f.396.0.2 <- NULL
bdE4$f.397.0.1 <- NULL
bdE4$f.397.0.2 <- NULL

# N. of incorrect matchings - Make it categorical #
bdE4$f.398.0.1 <- as.factor(bdE4$f.398.0.1)
bdE4$f.398.0.2 <- as.factor(bdE4$f.398.0.2)

# N. of incorrect matchings - exclude because inverse of correct #
bdE4$f.399.0.1  <- NULL
bdE4$f.399.0.2  <- NULL

# EXCLUSION SEVERAL VARIABLES IN SNAP GAME #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("401","402","403","404")]

# Own or rent accomodation lived in #
levels(bdE4$f.680.0.0) <- c(levels(bdE4$f.680.0.0),"No home")
bdE4$f.680.0.0[bdE4$f.670.0.0=="Sheltered accommodation" | bdE4$f.670.0.0=="Care home"] <- "No home"

# Length of time at current address #
bdE4$f.699.0.0 <- ifelse(bdE4$f.699.0.0==-10, 0, ifelse(bdE4$f.699.0.0==-1,NA,ifelse(bdE4$f.699.0.0==-3,NA,bdE4$f.699.0.0)))

# Number in the household #
bdE4$f.709.0.0b <- bdE4$f.709.0.0
bdE4$f.709.0.0 <- ifelse(bdE4$f.709.0.0==-1, NA, ifelse(bdE4$f.709.0.0==-3,NA,bdE4$f.709.0.0))
bdE4$f.709.0.0 <- as.factor(ifelse(bdE4$f.670.0.0=="Sheltered accommodation" | bdE4$f.670.0.0=="Care home","No home",as.character(cut2(bdE4$f.709.0.0))))

# N. of vehicles in household #
levels(bdE4$f.728.0.0) <- c(levels(bdE4$f.728.0.0),"No home")
bdE4$f.728.0.0[bdE4$f.670.0.0=="Sheltered accommodation" | bdE4$f.670.0.0=="Care home"] <- "No home"

# Average total household income before tax #
levels(bdE4$f.738.0.0) <- c(levels(bdE4$f.738.0.0),"No home")
bdE4$f.738.0.0[bdE4$f.670.0.0=="Sheltered accommodation" | bdE4$f.670.0.0=="Care home"] <- "No home"

# Time employed in main current job # 
bdE4$f.757.0.0 <- ifelse(bdE4$f.757.0.0==-10, 0, ifelse(bdE4$f.757.0.0==-1,NA,ifelse(bdE4$f.757.0.0==-3,NA,bdE4$f.757.0.0)))
bdE4$f.757.0.0 <- as.factor(ifelse(bdE4$f.6142.0.l2=="In paid employment or self-employed",as.character(cut2(bdE4$f.757.0.0)),"Not working"))

# Length of workign week for main job #
bdE4$f.767.0.0 <- ifelse(bdE4$f.767.0.0==-1, NA, ifelse(bdE4$f.767.0.0==-3,NA,bdE4$f.767.0.0))
bdE4$f.767.0.0 <- as.factor(ifelse(bdE4$f.6142.0.l2=="In paid employment or self-employed",as.character(cut2(bdE4$f.767.0.0)),"Not working"))
	
# Frequency for travelling from home to job workplace #
bdE4$f.777.0.0 <- ifelse(bdE4$f.777.0.0==-10, 0, ifelse(bdE4$f.777.0.0==-1,NA,ifelse(bdE4$f.777.0.0==-3,NA,bdE4$f.777.0.0)))
bdE4$f.777.0.0 <- as.factor(ifelse(bdE4$f.6142.0.l2=="In paid employment or self-employed",as.character(cut2(bdE4$f.777.0.0)),"Not working"))

# Distance between home and job place #
bdE4$f.796.0.0 <- ifelse(bdE4$f.796.0.0==-10, 0, ifelse(bdE4$f.796.0.0==-1,NA,ifelse(bdE4$f.796.0.0==-3,NA,bdE4$f.796.0.0)))
bdE4$f.796.0.0 <- as.factor(ifelse(bdE4$f.6142.0.l2=="In paid employment or self-employed",as.character(cut2(bdE4$f.796.0.0)),"Not working"))

# Job involve mainly walking or standing #
levels(bdE4$f.806.0.0) <- c(levels(bdE4$f.806.0.0),"Not working")
bdE4$f.806.0.0[bdE4$f.6142.0.l2==0] <- "Not working"

# Job involves heavy manual or physical work #
levels(bdE4$f.816.0.0) <- c(levels(bdE4$f.816.0.0),"Not working")
bdE4$f.816.0.0[bdE4$f.6142.0.l2==0] <- "Not working"

# Job involves shift work #
levels(bdE4$f.826.0.0) <- c(levels(bdE4$f.826.0.0),"Not working")
bdE4$f.826.0.0[bdE4$f.6142.0.l2==0] <- "Not working"

# Age completed full education #
bdE4$f.845.0.0 <-  ifelse(bdE4$f.845.0.0==-2, 0, ifelse(bdE4$f.845.0.0==-1,NA,ifelse(bdE4$f.845.0.0==-3,NA,bdE4$f.845.0.0)))
bdE4$f.845.0.0 <- as.factor(ifelse(bdE4$f.6138.0.l2==1,"College or University",as.character(cut2(bdE4$f.845.0.0))))

# N. of days/week walked 10+ minutes #
bdE4$f.864.0.0 <-  ifelse(bdE4$f.864.0.0==-1, NA,ifelse(bdE4$f.864.0.0==-3,NA,bdE4$f.864.0.0))
bdE4$f.864.0.0 <- as.factor(bdE4$f.864.0.0)

# Duration of walks #
bdE4$f.874.0.0 <-  ifelse(bdE4$f.874.0.0==-1, NA, ifelse(bdE4$f.874.0.0==-3,NA,bdE4$f.874.0.0))
bdE4$f.874.0.0 <- as.factor(ifelse(bdE4$f.864.0.0!=-2 & bdE4$f.864.0.0!=0,as.character(cut2(bdE4$f.874.0.0)),"Not walking"))

print("second")

# Days of moderate physical activity #
bdE4$f.884.0.0 <-  ifelse(bdE4$f.884.0.0==-1, NA, ifelse(bdE4$f.884.0.0==-3,NA,bdE4$f.884.0.0))

# Duration of moderate physical activity #
bdE4$f.894.0.0 <-  ifelse(bdE4$f.894.0.0==-1, NA, ifelse(bdE4$f.894.0.0==-3,NA,bdE4$f.894.0.0))
bdE4$f.894.0.0 <- as.factor(ifelse(bdE4$f.884.0.0 >0,as.character(cut2(bdE4$f.894.0.0)),"Not walking"))

# Days of vigorous physical activity #
bdE4$f.904.0.0 <-  ifelse(bdE4$f.904.0.0==-1, NA, ifelse(bdE4$f.904.0.0==-3,NA,bdE4$f.904.0.0))

# Duration of vigorous physical activity #
bdE4$f.914.0.0 <-  ifelse(bdE4$f.914.0.0==-1, NA, ifelse(bdE4$f.914.0.0==-3,NA,bdE4$f.914.0.0))
bdE4$f.914.0.0 <- as.factor(ifelse(bdE4$f.904.0.0>0,as.character(cut2(bdE4$f.914.0.0)),"Not walking"))

# Usual walking pace #
levels(bdE4$f.924.0.0) <- c(levels(bdE4$f.924.0.0),"Not walking")
bdE4$f.924.0.0[bdE4$f.864.0.0==-2] <- "Not walking"

# Frequency stairs climbing #
levels(bdE4$f.943.0.0) <- c(levels(bdE4$f.943.0.0),"Not walking")
bdE4$f.943.0.0[bdE4$f.864.0.0==-2] <- "Not walking"

# Frequency walking for pleasure #
levels(bdE4$f.971.0.0) <- c(levels(bdE4$f.971.0.0),"Not walking for pleasure")
bdE4$f.971.0.0[bdE4$f.6164.0.l2==0] <- "Not walking for pleasure"

# Duration walking for pleasure #
levels(bdE4$f.981.0.0) <- c(levels(bdE4$f.981.0.0),"Not walking for pleasure")
bdE4$f.981.0.0[bdE4$f.6164.0.l2==0] <- "Not walking for pleasure"

# Frequency of strenuous sports #
levels(bdE4$f.991.0.0) <- c(levels(bdE4$f.991.0.0),"Not doing strenuous sports")
bdE4$f.991.0.0[bdE4$f.6164.0.l4 == 0] <- "Not doing strenuous sports"

# Duration of strenuous sports #
levels(bdE4$f.1001.0.0) <- c(levels(bdE4$f.1001.0.0),"Not doing strenuous sports")
bdE4$f.1001.0.0[bdE4$f.6164.0.l4 == 0] <- "Not doing strenuous sports"

# Frequency of light DIY #
levels(bdE4$f.1011.0.0) <- c(levels(bdE4$f.1011.0.0),"Not light DIY")
bdE4$f.1011.0.0[bdE4$f.6164.0.l5 == 0] <- "Not light DIY"

# Duration of light DIY #
levels(bdE4$f.1021.0.0) <- c(levels(bdE4$f.1021.0.0),"Not light DIY")
bdE4$f.1021.0.0[bdE4$f.6164.0.l5 == 0 ] <- "Not light DIY"

# Time spend outdoors in the summer #
bdE4$f.1050.0.0 <- ifelse(bdE4$f.1050.0.0==-10, 0, ifelse(bdE4$f.1050.0.0==-1,NA,ifelse(bdE4$f.1050.0.0==-3,NA,bdE4$f.1050.0.0)))

# Time spend outdoors in the winter #
bdE4$f.1060.0.0 <- ifelse(bdE4$f.1060.0.0==-10, 0, ifelse(bdE4$f.1060.0.0==-1,NA,ifelse(bdE4$f.1060.0.0==-3,NA,bdE4$f.1060.0.0)))

# Time spend watching television #
bdE4$f.1070.0.0 <- ifelse(bdE4$f.1070.0.0==-10, 0, ifelse(bdE4$f.1070.0.0==-1,NA,ifelse(bdE4$f.1070.0.0==-3,NA,bdE4$f.1070.0.0)))

# Time spend using computer #
bdE4$f.1080.0.0 <- ifelse(bdE4$f.1080.0.0==-10, 0, ifelse(bdE4$f.1080.0.0==-1,NA,ifelse(bdE4$f.1080.0.0==-3,NA,bdE4$f.1080.0.0)))

# Time spend driving #
bdE4$f.1090.0.0 <- ifelse(bdE4$f.1090.0.0==-10, 0, ifelse(bdE4$f.1090.0.0==-1,NA,ifelse(bdE4$f.1090.0.0==-3,NA,bdE4$f.1090.0.0)))

# Weekly usage of mobile phone #
levels(bdE4$f.1120.0.0) <- c(levels(bdE4$f.1120.0.0),"No use phone")
bdE4$f.1120.0.0[bdE4$f.1110.0.0=="Never used mobile phone at least once per week"] <- "No use phone"

# Hand-free device #
levels(bdE4$f.1130.0.0) <- c(levels(bdE4$f.1130.0.0),"No use phone")
bdE4$f.1130.0.0[bdE4$f.1110.0.0=="Never used mobile phone at least once per week"] <- "No use phone"

# Difference compared to two years ago #
levels(bdE4$f.1140.0.0) <- c(levels(bdE4$f.1140.0.0),"No use phone")
bdE4$f.1140.0.0[bdE4$f.1110.0.0=="Never used mobile phone at least once per week"] <- "No use phone"

# Side head #
levels(bdE4$f.1150.0.0) <- c(levels(bdE4$f.1150.0.0),"No use phone")
bdE4$f.1150.0.0[bdE4$f.1110.0.0=="Never used mobile phone at least once per week"] <- "No use phone"

# Sleep duration #
bdE4$f.1160.0.0 <-  ifelse(bdE4$f.1160.0.0==-1, NA, ifelse(bdE4$f.1160.0.0==-3,NA,bdE4$f.1160.0.0))

# Past tobacco smoking #
levels(bdE4$f.1249.0.0) <- c(levels(bdE4$f.1249.0.0),"I'm a smoker")
bdE4$f.1249.0.0[bdE4$f.1239.0.0=="Yes, on most or all days"] <- "I'm a smoker"

# Smoker in household #
levels(bdE4$f.1259.0.0) <- c(levels(bdE4$f.1259.0.0),"I'm a smoker")
bdE4$f.1259.0.0[bdE4$f.1239.0.0=="Yes, on most or all days"] <- "I'm a smoker"

# Exposure to tabacco smoke at home - make two categories 0, >0 #
bdE4$f.1269.0.0 <-  as.factor(ifelse(bdE4$f.1269.0.0==-1, NA, ifelse(bdE4$f.1269.0.0==-3,NA,ifelse(bdE4$f.1269.0.0==0,0,">=1"))))
levels(bdE4$f.1269.0.0) <- c(levels(bdE4$f.1269.0.0),"I'm a smoker")
bdE4$f.1269.0.0[bdE4$f.1239.0.0=="Yes, on most or all days"] <- "I'm a smoker"

# Exposure to tabacco outside - make two categories 0, >0 #
bdE4$f.1279.0.0 <-  as.factor(ifelse(bdE4$f.1279.0.0==-1, NA, ifelse(bdE4$f.1279.0.0==-3,NA,ifelse(bdE4$f.1279.0.0==0,0,1))))
levels(bdE4$f.1279.0.0) <- c(levels(bdE4$f.1279.0.0),"I'm a smoker")
bdE4$f.1279.0.0[bdE4$f.1239.0.0=="Yes, on most or all days"] <- "I'm a smoker"

# Cooked vegetable intake #
bdE4$f.1289.0.0 <-  ifelse(bdE4$f.1289.0.0==-10, 0, ifelse(bdE4$f.1289.0.0==-1,NA,ifelse(bdE4$f.1289.0.0==-3,NA,bdE4$f.1289.0.0)))

# Raw vegetable intake #
bdE4$f.1299.0.0 <-  ifelse(bdE4$f.1299.0.0==-10, 0, ifelse(bdE4$f.1299.0.0==-1,NA,ifelse(bdE4$f.1299.0.0==-3,NA,bdE4$f.1299.0.0)))

# Fresh fruit intake #
bdE4$f.1309.0.0 <-  ifelse(bdE4$f.1309.0.0==-10, 0, ifelse(bdE4$f.1309.0.0==-1,NA,ifelse(bdE4$f.1309.0.0==-3,NA,bdE4$f.1309.0.0)))

# Dried fruit intake #
bdE4$f.1319.0.0 <-  ifelse(bdE4$f.1319.0.0==-10, 0, ifelse(bdE4$f.1319.0.0==-1,NA,ifelse(bdE4$f.1319.0.0==-3,NA,bdE4$f.1319.0.0)))

# Cheese intake #
levels(bdE4$f.1408.0.0) <- c(levels(bdE4$f.1408.0.0),"Not eat cheese")
bdE4$f.1408.0.0[bdE4$f.6144.0.l2==1] <- "Not eat cheese"

# Bread intake #
bdE4$f.1438.0.0b <- bdE4$f.1438.0.0
bdE4$f.1438.0.0 <-  ifelse(bdE4$f.1438.0.0==-10, 0, ifelse(bdE4$f.1438.0.0==-1,NA,ifelse(bdE4$f.1438.0.0==-3,NA,bdE4$f.1438.0.0)))

# Bread type #
levels(bdE4$f.1448.0.0) <- c(levels(bdE4$f.1448.0.0),"Not eat bread")
bdE4$f.1448.0.0[bdE4$f.1438.0.0b<1] <- "Not eat bread"

# Cereal intake #
bdE4$f.1458.0.0b <- bdE4$f.1458.0.0
bdE4$f.1458.0.0 <-  ifelse(bdE4$f.1458.0.0==-10, 0, ifelse(bdE4$f.1458.0.0==-1,NA,ifelse(bdE4$f.1458.0.0==-3,NA,bdE4$f.1458.0.0)))

# Cereal type #
levels(bdE4$f.1468.0.0) <- c(levels(bdE4$f.1468.0.0),"Not eat cereal")
bdE4$f.1468.0.0[bdE4$f.1458.0.0b<1] <- "Not eat cereal"

# Tea intake #
bdE4$f.1488.0.0 <-  ifelse(bdE4$f.1488.0.0==-10, 0, ifelse(bdE4$f.1488.0.0==-1,NA,ifelse(bdE4$f.1488.0.0==-3,NA,bdE4$f.1488.0.0)))

# Coffee intake #
bdE4$f.1498.0.0 <-  ifelse(bdE4$f.1498.0.0==-10, 0, ifelse(bdE4$f.1498.0.0==-1,NA,ifelse(bdE4$f.1498.0.0==-3,NA,bdE4$f.1498.0.0)))
bdE4$f.1498.0.0b <- bdE4$f.1498.0.0

# Coffee type #
levels(bdE4$f.1508.0.0) <- c(levels(bdE4$f.1508.0.0),"Not drinking coffee")
bdE4$f.1508.0.0[bdE4$f.1498.0.0b<2] <- "Not drinking coffee"

# Water intake #
bdE4$f.1528.0.0 <-  ifelse(bdE4$f.1528.0.0==-10, 0, ifelse(bdE4$f.1528.0.0==-1,NA,ifelse(bdE4$f.1528.0.0==-3,NA,bdE4$f.1528.0.0)))

# Weekly red wine intake #
bdE4$f.1568.0.0 <-  ifelse(bdE4$f.1568.0.0==-1, NA, ifelse(bdE4$f.1568.0.0==-3,NA,bdE4$f.1568.0.0))
bdE4$f.1568.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Daily or almost daily" | bdE4$f.1558.0.0=="Three or four times a week" | bdE4$f.1558.0.0=="Once or twice a week",as.character(cut2(bdE4$f.1568.0.0)),"Not drinking / Not drinking often"))

# Weekly champagne + white wine intake #
bdE4$f.1578.0.0 <-  ifelse(bdE4$f.1578.0.0==-1, NA, ifelse(bdE4$f.1578.0.0==-3,NA,bdE4$f.1578.0.0))
bdE4$f.1578.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Daily or almost daily" | bdE4$f.1558.0.0=="Three or four times a week" | bdE4$f.1558.0.0=="Once or twice a week",as.character(cut2(bdE4$f.1578.0.0)),"Not drinking / Not drinking often"))

# Weekly beer + cider intake #
bdE4$f.1588.0.0 <-  ifelse(bdE4$f.1588.0.0==-1, NA, ifelse(bdE4$f.1588.0.0==-3,NA,bdE4$f.1588.0.0))
bdE4$f.1588.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Daily or almost daily" | bdE4$f.1558.0.0=="Three or four times a week" | bdE4$f.1558.0.0=="Once or twice a week",as.character(cut2(bdE4$f.1588.0.0)),"Not drinking / Not drinking often"))

# Weekly spirits + cider intake - Make two categories 0, >0 #
bdE4$f.1598.0.0 <-  as.factor(ifelse(bdE4$f.1598.0.0==-1, NA, ifelse(bdE4$f.1598.0.0==-3,NA,ifelse(bdE4$f.1598.0.0==0,0,">=1"))))
levels(bdE4$f.1598.0.0) <- c(levels(bdE4$f.1598.0.0),"Not drinking / Not drinking often")
bdE4$f.1598.0.0[bdE4$f.1558.0.0!="Daily or almost daily" & bdE4$f.1558.0.0!="Three or four times a week" & bdE4$f.1558.0.0!="Once or twice a week"] <- "Not drinking / Not drinking often"

# Weekly fortified wine intake - Make two categories 0, >0 #
bdE4$f.1608.0.0 <-  as.factor(ifelse(bdE4$f.1608.0.0==-1, NA, ifelse(bdE4$f.1608.0.0==-3,NA,ifelse(bdE4$f.1608.0.0==0,0,">=1"))))
levels(bdE4$f.1608.0.0) <- c(levels(bdE4$f.1608.0.0),"Not drinking / Not drinking often")
bdE4$f.1608.0.0[bdE4$f.1558.0.0!="Daily or almost daily" & bdE4$f.1558.0.0!="Three or four times a week" & bdE4$f.1558.0.0!="Once or twice a week"] <- "Not drinking / Not drinking often"

# Alchool during meals #
levels(bdE4$f.1618.0.0) <- c(levels(bdE4$f.1618.0.0),"Not drinking")
bdE4$f.1618.0.0[bdE4$f.1558.0.0=="Never"] <- "Not drinking"

# Alchool compared 10-years #
levels(bdE4$f.1628.0.0) <- c(levels(bdE4$f.1628.0.0),"Not drinking")
bdE4$f.1628.0.0[bdE4$f.1558.0.0=="Never"] <- "Not drinking"

# Child subburn occasion #
bdE4$f.1737.0.0 <-  ifelse(bdE4$f.1737.0.0==-1, NA, ifelse(bdE4$f.1737.0.0==-3,NA,bdE4$f.1737.0.0))

# Part multiple birth #
levels(bdE4$f.1777.0.0) <- c(levels(bdE4$f.1777.0.0),"Adopdted")
bdE4$f.1777.0.0[bdE4$f.1767.0.0=="Yes"] <- "Adopdted"

# Maternal smoking #
levels(bdE4$f.1787.0.0) <- c(levels(bdE4$f.1787.0.0),"Adopdted")
bdE4$f.1787.0.0[bdE4$f.1767.0.0=="Yes"] <- "Adopdted"

# Father still alive #
levels(bdE4$f.1797.0.0) <- c(levels(bdE4$f.1797.0.0),"Adopdted")
bdE4$f.1797.0.0[bdE4$f.1767.0.0=="Yes"] <- "Adopdted"

# Father's age death #
bdE4$f.1807.0.0 <-  ifelse(bdE4$f.1807.0.0==-1, NA, ifelse(bdE4$f.1807.0.0==-3,NA,bdE4$f.1807.0.0))
bdE4$f.1807.0.0 <- as.factor(ifelse(bdE4$f.1797.0.0=="No",as.character(cut2(bdE4$f.1807.0.0)),"Still alive"))

# Mother still alive #
levels(bdE4$f.1835.0.0) <- c(levels(bdE4$f.1835.0.0),"Adopdted")
bdE4$f.1835.0.0[bdE4$f.1767.0.0=="Yes"] <- "Adopdted"

# Mother's age #
bdE4$f.1845.0.0 <-  ifelse(bdE4$f.1845.0.0==-1, NA, ifelse(bdE4$f.1845.0.0==-3,NA,bdE4$f.1845.0.0))
bdE4$f.1845.0.0 <- as.factor(ifelse(bdE4$f.1835.0.0=="Yes",as.character(cut2(bdE4$f.1845.0.0)),"Not alive"))

# N. of full brothers #
bdE4$f.1873.0.0b <- bdE4$f.1873.0.0
bdE4$f.1873.0.0 <-  ifelse(bdE4$f.1873.0.0==-1, NA, ifelse(bdE4$f.1873.0.0==-3,NA,bdE4$f.1873.0.0))
bdE4$f.1873.0.0 <- as.factor(ifelse(bdE4$f.1767.0.0=="Yes","Adopted",as.character(cut2(bdE4$f.1873.0.0))))

# N. of full sisters #
bdE4$f.1883.0.0b <- bdE4$f.1883.0.0
bdE4$f.1883.0.0 <-  ifelse(bdE4$f.1883.0.0==-1, NA, ifelse(bdE4$f.1883.0.0==-3,NA,bdE4$f.1883.0.0))
bdE4$f.1883.0.0 <- as.factor(ifelse(bdE4$f.1767.0.0=="Yes","Adopted",as.character(cut2(bdE4$f.1883.0.0))))

# Age first sexual intercurse #
bdE4$f.2139.0.0 <-  ifelse(bdE4$f.2139.0.0==-1, NA,ifelse(bdE4$f.2139.0.0==-3,NA,bdE4$f.2139.0.0))
bdE4$f.2139.0.0 <- as.factor(ifelse(bdE4$f.2139.0.0==-2,"No sex",as.character(cut2(bdE4$f.2139.0.0))))

# Lifetime number of sexual partners #
bdE4$f.2149.0.0 <-  ifelse(bdE4$f.2149.0.0==-1, NA, ifelse(bdE4$f.2149.0.0==-2,0,ifelse(bdE4$f.2149.0.0==-3,NA,bdE4$f.2149.0.0)))
bdE4$f.2149.0.0 <- as.factor(ifelse(bdE4$f.2139.0.0=="No sex","No sex",as.character(cut2(bdE4$f.2149.0.0))))

# Same sex intercurse #
levels(bdE4$f.2159.0.0) <- c(levels(bdE4$f.2159.0.0),"No sex")
bdE4$f.2159.0.0[bdE4$f.2139.0.0=="No sex"] <- "No sex"

# Age started to wear glasses or contact lenses
bdE4$f.2217.0.0 <-  ifelse(bdE4$f.2217.0.0==-1, NA, ifelse(bdE4$f.2217.0.0==-3,NA,bdE4$f.2217.0.0))
bdE4$f.2217.0.0 <- as.factor(ifelse(bdE4$f.2207.0.0=="No","No glasses",as.character(cut2(bdE4$f.2217.0.0))))

# Hearing difficulty/problems with background noise #
levels(bdE4$f.2257.0.0) <- c(levels(bdE4$f.2257.0.0),"Deaf")
bdE4$f.2257.0.0[bdE4$f.2247.0.0=="I am completely deaf"] <- "Deaf"

# Frequency of solarium/sunlamp use - Make three categories, 0,1-5,>5 #
bdE4$f.2277.0.0 <- as.factor(ifelse(bdE4$f.2277.0.0==-10, 0, ifelse(bdE4$f.2277.0.0==-1,NA,ifelse(bdE4$f.2277.0.0==-3,NA,ifelse(bdE4$f.2277.0.0==0,0,ifelse(bdE4$f.2277.0.0<10,"1-10",">10"))))))

# Most recent bowel cancer screening #
bdE4$f.2355.0.0 <- ifelse(bdE4$f.2355.0.0==-10, 0, ifelse(bdE4$f.2355.0.0==-1,NA,ifelse(bdE4$f.2355.0.0==-3,NA,bdE4$f.2355.0.0)))
bdE4$f.2355.0.0 <- as.factor(ifelse(bdE4$f.2345.0.0=="Yes",as.character(cut2(bdE4$f.2355.0.0)),"No screening"))


# N. of children fathered #
bdE4$f.2405.0.0 <-  ifelse(bdE4$f.2405.0.0==-1, NA, ifelse(bdE4$f.2405.0.0==-3,NA,bdE4$f.2405.0.0))
bdE4$f.2405.0.0 <- cut2(bdE4$f.2405.0.0)


# Frequency of heavy DIY #
levels(bdE4$f.2624.0.0) <- c(levels(bdE4$f.2624.0.0),"Not heavy DIY")
bdE4$f.2624.0.0[bdE4$f.6164.0.l6==0] <- "Not heavy DIY"

# Duration of heavy DIY #
levels(bdE4$f.2634.0.0) <- c(levels(bdE4$f.2634.0.0),"Not heavy DIY")
bdE4$f.2634.0.0[bdE4$f.6164.0.l6==0] <- "Not heavy DIY"

# Light smokers #
levels(bdE4$f.2644.0.0) <- c(levels(bdE4$f.2644.0.0),"Not smoking / No light smokers")
bdE4$f.2644.0.0[bdE4$f.1249.0.0!="Smoked occasionally" & bdE4$f.1249.0.0!="Just tried once or twice"] <- "Not smoking / No light smokers"

# Non-butter spread type #
levels(bdE4$f.2654.0.0) <- c(levels(bdE4$f.2654.0.0),"No butter")
bdE4$f.2654.0.0[bdE4$f.1428.0.0!="Other type of spread/margarine" & bdE4$f.1428.0.0!="Do not know"] <- "No butter"

# Reason to reduce alchool #
levels(bdE4$f.2664.0.0) <- c(levels(bdE4$f.2664.0.0),"Not reduced or never drinking")
bdE4$f.2664.0.0[bdE4$f.1558.0.0=="Never" | bdE4$f.1628.0.0!="Less nowadays" ] <- "Not reduced or never drinking"


# Years since last breast cancer screening #
bdE4$f.2684.0.0 <- ifelse(bdE4$f.2684.0.0==-10, 0, ifelse(bdE4$f.2684.0.0==-1,NA,ifelse(bdE4$f.2684.0.0==-3,NA,bdE4$f.2684.0.0)))
bdE4$f.2684.0.0 <- as.factor(ifelse(bdE4$f.2674.0.0=="No" ,"No screening",as.character(cut2(bdE4$f.2684.0.0))))

# Years since last cervical smear test #
bdE4$f.2704.0.0 <- ifelse(bdE4$f.2704.0.0==-10, 0, ifelse(bdE4$f.2704.0.0==-1,NA,ifelse(bdE4$f.2704.0.0==-3,NA,bdE4$f.2704.0.0)))
bdE4$f.2704.0.0 <- as.factor(ifelse(bdE4$f.2694.0.0=="No","No screening",as.character(cut2(bdE4$f.2704.0.0))))

# Age when period started #
bdE4$f.2714.0.0 <-  ifelse(bdE4$f.2714.0.0==-1, NA, ifelse(bdE4$f.2714.0.0==-3,NA,bdE4$f.2714.0.0))
bdE4$f.2714.0.0 <- cut2(bdE4$f.2714.0.0)

# N. of life births - b version is used later #
bdE4$f.2734.0.0 <-  ifelse(bdE4$f.2734.0.0==-1, NA, ifelse(bdE4$f.2734.0.0==-3,NA,bdE4$f.2734.0.0))
bdE4$f.2734.0.0b <- bdE4$f.2734.0.0
bdE4$f.2734.0.0 <- cut2(bdE4$f.2734.0.0)

# Birth weight of first child #
bdE4$f.2744.0.0 <-  ifelse(bdE4$f.2744.0.0==-1, NA,ifelse(bdE4$f.2744.0.0==-3,NA,bdE4$f.2744.0.0))
bdE4$f.2744.0.0 <- as.factor(ifelse(bdE4$f.2734.0.0b < 1 ,"No child",ifelse(bdE4$f.2744.0.0==-2,"Twins",as.character(cut2(bdE4$f.2744.0.0)))))

# Age of first birth #
bdE4$f.2754.0.0 <-  ifelse(bdE4$f.2754.0.0==-4, NA, ifelse(bdE4$f.2754.0.0==-3,NA,bdE4$f.2754.0.0))
bdE4$f.2754.0.0 <- as.factor(ifelse(bdE4$f.2734.0.0b < 2 ,"Less than 2 children",as.character(cut2(bdE4$f.2754.0.0))))

# Age of last birth #
bdE4$f.2764.0.0 <-  ifelse(bdE4$f.2764.0.0==-4, NA, ifelse(bdE4$f.2764.0.0==-3,NA,bdE4$f.2764.0.0))
bdE4$f.2764.0.0 <- as.factor(ifelse(bdE4$f.2734.0.0b < 2,"Less than 2 children",as.character(cut2(bdE4$f.2764.0.0))))

# Age started oral contraceptive pill #
bdE4$f.2794.0.0 <-  ifelse(bdE4$f.2794.0.0== -1, NA, ifelse(bdE4$f.2794.0.0== -3,NA,bdE4$f.2794.0.0))
bdE4$f.2794.0.0 <- as.factor(ifelse(bdE4$f.2784.0.0 == "No","No pill",as.character(cut2(bdE4$f.2794.0.0))))

# Age last used oral contraceptive pill #
bdE4$f.2804.0.0b <- bdE4$f.2804.0.0
bdE4$f.2804.0.0 <-  ifelse(bdE4$f.2804.0.0 == -1, NA, ifelse(bdE4$f.2804.0.0 == -3,NA,bdE4$f.2804.0.0))
bdE4$f.2804.0.0 <- as.factor(ifelse(bdE4$f.2784.0.0 == "No","No pill",ifelse(bdE4$f.2804.0.0b == -11,"Still taking",as.character(cut2(bdE4$f.2804.0.0)))))


# Age at hysterectomy #
bdE4$f.2824.0.0 <-  ifelse(bdE4$f.2824.0.0==-1, NA, ifelse(bdE4$f.2824.0.0==-3,NA,bdE4$f.2824.0.0))
bdE4$f.2824.0.0 <- as.factor(ifelse((bdE4$f.2724.0.0 != "Not sure - had a hysterectomy" & bdE4$f.3591.0.0!="Yes"),"No hysterectomy",as.character(cut2(bdE4$f.2824.0.0))))


# Age smoking in former smokers #
bdE4$f.2867.0.0 <-  ifelse(bdE4$f.2867.0.0==-1, NA, ifelse(bdE4$f.2867.0.0==-3,NA,bdE4$f.2867.0.0))
bdE4$f.2867.0.0 <- as.factor(ifelse(bdE4$f.1249.0.0 == "Smoked on most or all days",as.character(cut2(bdE4$f.2867.0.0)),"No past smokers"))

# Type of tobacco previously smoked #
levels(bdE4$f.2877.0.0) <- c(levels(bdE4$f.2877.0.0),"No past smokers")
bdE4$f.2877.0.0[bdE4$f.1249.0.0 != "Smoked on most or all days"] <- "No past smokers"

# N. cig previsouly smoked daily #
bdE4$f.2887.0.0 <-  ifelse(bdE4$f.2887.0.0==-10, 0, ifelse(bdE4$f.2887.0.0==-1,NA,bdE4$f.2887.0.0))
bdE4$f.2887.0.0 <- as.factor(ifelse(bdE4$f.1249.0.0 == "Smoked on most or all days" & (bdE4$f.2877.0.0=="Hand-rolled cigarettes" | bdE4$f.2877.0.0=="Manufactured cigarettes"),as.character(cut2(bdE4$f.2887.0.0)),"No past smokers or no cigarettes"))

# Age stop smoking #
bdE4$f.2897.0.0 <-  ifelse(bdE4$f.2897.0.0==-1, NA, ifelse(bdE4$f.2897.0.0==-3,NA,bdE4$f.2897.0.0))
bdE4$f.2897.0.0 <- as.factor(ifelse(bdE4$f.1249.0.0 == "Smoked on most or all days" ,as.character(cut2(bdE4$f.2897.0.0)),"No past smokers"))

# Ever stopped smoking for 6+ months #
levels(bdE4$f.2907.0.0) <- c(levels(bdE4$f.2907.0.0),"No past smokers")
bdE4$f.2907.0.0[bdE4$f.1249.0.0 != "Smoked on most or all days"] <- "No past smokers"

# N. of unsuccessful stop-smoking attempts #
bdE4$f.2926.0.0 <-  ifelse(bdE4$f.2926.0.0==-1, NA, ifelse(bdE4$f.2926.0.0==-3,NA,bdE4$f.2926.0.0))
bdE4$f.2926.0.0 <- as.factor(ifelse(bdE4$f.1249.0.0 == "Smoked on most or all days" ,as.character(cut2(bdE4$f.2926.0.0)),"No past smokers"))

# Likelihood of resuming smoking #
levels(bdE4$f.2936.0.0) <- c(levels(bdE4$f.2936.0.0),"No past smokers")
bdE4$f.2936.0.0[bdE4$f.1249.0.0 != "Smoked on most or all days"] <- "No past smokers"

# Father age #
bdE4$f.2946.0.0 <-  ifelse(bdE4$f.2946.0.0==-1, NA, ifelse(bdE4$f.2946.0.0==-3,NA,bdE4$f.2946.0.0))
bdE4$f.2946.0.0 <- as.factor(ifelse(bdE4$f.1797.0.0 == "Yes" ,as.character(cut2(bdE4$f.2946.0.0)),"Father dead"))

# General pain more 3 months #
levels(bdE4$f.2956.0.0) <- c(levels(bdE4$f.2956.0.0),"No pain all over the body")
bdE4$f.2956.0.0[bdE4$f.6159.0.l9 ==0] <- "No pain all over the body"

# Age high blood pressure diagnosed #
bdE4$f.2966.0.0 <-  ifelse(bdE4$f.2966.0.0==-1, NA, ifelse(bdE4$f.2966.0.0==-3,NA,bdE4$f.2966.0.0))
bdE4$f.2966.0.0 <- as.factor(ifelse(bdE4$f.6150.0.l5 == "High blood pressure" ,as.character(cut2(bdE4$f.2966.0.0)),"No High blood pressure"))

# Age diabetes diagnosed #
bdE4$f.2976.0.0 <-  ifelse(bdE4$f.2976.0.0==-1, NA, ifelse(bdE4$f.2976.0.0==-3,NA,bdE4$f.2976.0.0))
bdE4$f.2976.0.0 <- as.factor(ifelse(bdE4$f.2443.0.0 == "Yes" ,as.character(cut2(bdE4$f.2976.0.0)),"No diabetes"))

print("third")

# Started insulin within one year from diabetes #
levels(bdE4$f.2986.0.0) <- c(levels(bdE4$f.2986.0.0),"No diabetes")
bdE4$f.2986.0.0[bdE4$f.2443.0.0 != "Yes"] <- "No diabetes"

# Fracture from simple fall  #
levels(bdE4$f.3005.0.0) <- c(levels(bdE4$f.3005.0.0),"No fractures")
bdE4$f.3005.0.0[bdE4$f.2463.0.0 != "Yes"] <- "No fractures"

# Results ranking #
bdE4$f.3059.0.0 <- NULL
bdE4$f.3059.0.1 <- NULL
bdE4$f.3059.0.2 <- NULL

# Forced vital capacity - keep only first measurement  #
bdE4$f.3062.0.1 <- NULL
bdE4$f.3062.0.2 <- NULL

# Forced expiratory volume - keep only first measurement #
bdE4$f.3063.0.1 <- NULL
bdE4$f.3063.0.2 <- NULL

# Peak expiratory flow - keep only first measurement #
bdE4$f.3064.0.1 <- NULL
bdE4$f.3064.0.2 <- NULL

# Ordering of blows #
bdE4$f.3065.0.0 <- NULL
bdE4$f.3065.0.1 <- NULL
bdE4$f.3065.0.2 <- NULL

# Foot measured for bone density 
bdE4$f.3081.0.0 <- NULL

# Fractured heel 
bdE4$f.3082.0.0 <- NULL

# Number of measurements made #
bdE4$f.3137.0.0 <- NULL

# Heel Broadband ultrasound attenuation, direct entry #
bdE4$f.3144.0.0 <- ifelse(is.na(bdE4$f.3144.0.0),bdE4$f.3085.0.0,bdE4$f.3144.0.0)
bdE4$f.3085.0.0 <- NULL

# Heel quantitative ultrasound index, direct entry #
bdE4$f.3147.0.0 <- ifelse(is.na(bdE4$f.3147.0.0),bdE4$f.3083.0.0,bdE4$f.3147.0.0)
bdE4$f.3083.0.0 <- NULL

# Heel bone mineral density (BMD) #
bdE4$f.3148.0.0 <- ifelse(is.na(bdE4$f.3148.0.0),bdE4$f.3084.0.0,bdE4$f.3148.0.0)
bdE4$f.3084.0.0 <- NULL

# Contra-indication for spirometry, too few #
bdE4$f.3088.0.0 <- NULL

# Neck/shoulder pain for 3+ months #
levels(bdE4$f.3404.0.0) <- c(levels(bdE4$f.3404.0.0),"No pain shoulder/neck")
bdE4$f.3404.0.0[bdE4$f.6159.0.l4 == 0] <- "No pain shoulder/neck"

# Hip pain for 3+ months #
levels(bdE4$f.3414.0.0) <- c(levels(bdE4$f.3414.0.0),"No pain hip")
bdE4$f.3414.0.0[bdE4$f.6159.0.l7 == 0] <- "No pain hip"

# Job involving night shift #
levels(bdE4$f.3426.0.0) <- c(levels(bdE4$f.3426.0.0),"No night shift or not working")
bdE4$f.3426.0.0[bdE4$f.6142.0.l2 == 0 |  bdE4$f.826.0.0 == "Never/rarely"] <- "No night shift or not working"

# Age started smoking in current smokers #
bdE4$f.3436.0.0 <-  ifelse(bdE4$f.3436.0.0==-1, NA, ifelse(bdE4$f.3436.0.0==-3,NA,bdE4$f.3436.0.0))
bdE4$f.3436.0.0 <- as.factor(ifelse(bdE4$f.1239.0.0 == "Yes, on most or all days" ,as.character(cut2(bdE4$f.3436.0.0)),"Not current smokers"))

# Type of currently smoked tobacco #
levels(bdE4$f.3446.0.0) <- c(levels(bdE4$f.3446.0.0),"Not current smokers")
bdE4$f.3446.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days"] <- "Not current smokers"

# N. of cigarettes currently smoked 
bdE4$f.3456.0.0 <-  ifelse(bdE4$f.3456.0.0==-10, 0, ifelse(bdE4$f.3456.0.0==-1,NA,ifelse(bdE4$f.3456.0.0==-3,NA,bdE4$f.3456.0.0)))
bdE4$f.3456.0.0 <- as.factor(ifelse(bdE4$f.3446.0.0 == "Hand-rolled cigarettes" | bdE4$f.3446.0.0=="Manufactured cigarettes",as.character(cut2(bdE4$f.3456.0.0)),"Not current smokers or no cigarettes"))

# Time from waking to first cigarette #
levels(bdE4$f.3466.0.0) <- c(levels(bdE4$f.3466.0.0),"Not current smokers or no cigarettes")
bdE4$f.3466.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days" | (bdE4$f.3446.0.0!="Manufactured cigarettes" & bdE4$f.3446.0.0!="Hand-rolled cigarettes")] <- "Not current smokers or no cigarettes"

# Difficulty not smoking for 1 day #
levels(bdE4$f.3476.0.0) <- c(levels(bdE4$f.3476.0.0),"Not current smokers or no cigarettes")
bdE4$f.3476.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days" | (bdE4$f.3446.0.0!="Manufactured cigarettes" & bdE4$f.3446.0.0!="Hand-rolled cigarettes")] <- "Not current smokers or no cigarettes"

# Ever tried to stop smoking #
levels(bdE4$f.3486.0.0) <- c(levels(bdE4$f.3486.0.0),"Not current smokers")
bdE4$f.3486.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days"] <- "Not current smokers"

# Wants to stop smoking #
levels(bdE4$f.3496.0.0) <- c(levels(bdE4$f.3496.0.0),"Not current smokers")
bdE4$f.3496.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days"] <- "Not current smokers"

# Smoking compared with 10 years previsouly #
levels(bdE4$f.3506.0.0) <- c(levels(bdE4$f.3506.0.0),"Not current smokers")
bdE4$f.3506.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days"] <- "Not current smokers"

# Mother age at death #
bdE4$f.3526.0.0 <- ifelse(bdE4$f.3526.0.0==-1, NA, ifelse(bdE4$f.3526.0.0==-3,NA,bdE4$f.3526.0.0))
bdE4$f.3526.0.0 <- as.factor(ifelse(bdE4$f.1835.0.0=="No",as.character(cut2(bdE4$f.3526.0.0)),"Still alive"))

# Age started hormone-replacment therapy #
bdE4$f.3536.0.0 <- ifelse(bdE4$f.3536.0.0==-1, NA, ifelse(bdE4$f.3536.0.0==-3,NA,bdE4$f.3536.0.0))
bdE4$f.3536.0.0 <- as.factor(ifelse(bdE4$f.2814.0.0=="Yes" ,as.character(cut2(bdE4$f.3536.0.0)),"No hormone therapy"))

# Age last used hormone-replacment therapy #
bdE4$f.3546.0.0 <- ifelse(bdE4$f.3546.0.0==-1, NA, ifelse(bdE4$f.3546.0.0==-3,NA,bdE4$f.3546.0.0))
bdE4$f.3546.0.0 <- as.factor(ifelse(bdE4$f.2814.0.0=="No", "No hormons",ifelse(bdE4$f.3546.0.0== -11,"Still taking", as.character(cut2(bdE4$f.3546.0.0)))))

# Back pain for 3+ months #
levels(bdE4$f.3571.0.0) <- c(levels(bdE4$f.3571.0.0),"No back pain")
bdE4$f.3571.0.0[bdE4$f.6159.0.l5 == 0] <- "No back pain"

# Age at menopause #
bdE4$f.3581.0.0 <- ifelse(bdE4$f.3581.0.0==-1, NA, ifelse(bdE4$f.3581.0.0==-3,NA,bdE4$f.3581.0.0))
bdE4$f.3581.0.0 <- as.factor(ifelse(bdE4$f.2724.0.0=="Yes" ,as.character(cut2(bdE4$f.3581.0.0)),"No menopause"))

# Ever had hysterectomy #
levels(bdE4$f.3591.0.0) <- c(levels(bdE4$f.3591.0.0),"Hysterectomy (in another question)")
bdE4$f.3591.0.0[bdE4$f.2724.0.0=="Not sure - had a hysterectomy"] <- "Hysterectomy (in another question)"

# Chest pain #
levels(bdE4$f.3606.0.0) <- c(levels(bdE4$f.3606.0.0),"No chest pain")
bdE4$f.3606.0.0[bdE4$f.2335.0.0 != "Yes"] <- "No chest pain"

# Chest pain ceases when still #
levels(bdE4$f.3616.0.0) <- c(levels(bdE4$f.3616.0.0),"No chest pain")
bdE4$f.3616.0.0[bdE4$f.2335.0.0 != "Yes" | (bdE4$f.3751.0.0!="Yes" & bdE4$f.3606.0.0!="Yes")] <- "No chest pain"

# Age at angina #
bdE4$f.3627.0.0 <- ifelse(bdE4$f.3627.0.0==-1, NA, ifelse(bdE4$f.3627.0.0==-3,NA,bdE4$f.3627.0.0))
bdE4$f.3627.0.0 <- as.factor(ifelse(bdE4$f.6150.0.l3=="Angina" ,as.character(cut2(bdE4$f.3627.0.0)),"No Angina"))

# Frequency of other exercises #
levels(bdE4$f.3637.0.0) <- c(levels(bdE4$f.3637.0.0),"No other exercise")
bdE4$f.3637.0.0[bdE4$f.6164.0.l3 == 0] <- "No other exercise"

# Duration of other exercise #
levels(bdE4$f.3647.0.0) <- c(levels(bdE4$f.3647.0.0),"No other exercise")
bdE4$f.3647.0.0[bdE4$f.6164.0.l3 == 0] <- "No other exercise"

# Year immigrated in UK #
bdE4$f.3659.0.0 <- ifelse(bdE4$f.3659.0.0==-1, NA, ifelse(bdE4$f.3659.0.0==-3,NA,bdE4$f.3659.0.0))
bdE4$f.3659.0.0 <- as.factor(ifelse(bdE4$f.1647.0.0=="Elsewhere" ,as.character(cut2(bdE4$f.3659.0.0)),"Born in UK"))

# Lifetime number of same-sex sexual partners #
bdE4$f.3669.0.0 <- ifelse(bdE4$f.3669.0.0==-1, NA, ifelse(bdE4$f.3669.0.0==-3,NA,bdE4$f.3669.0.0))
bdE4$f.3669.0.0 <- as.factor(ifelse(bdE4$f.2159.0.0=="Yes" ,as.character(cut2(bdE4$f.3669.0.0)),"No same-sex partner"))

# Age when last ate meat - exclude because no clear inclusion definition #
bdE4$f.3680.0.0 <- NULL

# Time since last menstrual period #
bdE4$f.3700.0.0 <- ifelse(bdE4$f.3700.0.0==-1, NA, ifelse(bdE4$f.3700.0.0==-3,NA,bdE4$f.3700.0.0))
bdE4$f.3700.0.0 <- as.factor(ifelse(bdE4$f.2724.0.0=="No",as.character(cut2(bdE4$f.3700.0.0)),"No menopause"))

# Length of menstrual period #
bdE4$f.3710.0.0b <- bdE4$f.3710.0.0
bdE4$f.3710.0.0 <- ifelse(bdE4$f.3710.0.0==-1, NA, ifelse(bdE4$f.3710.0.0==-3,NA,bdE4$f.3710.0.0))
bdE4$f.3710.0.0 <- as.factor(ifelse(bdE4$f.2724.0.0=="Yes","Not menopause",ifelse(bdE4$f.3710.0.0b==-6,"Irregular",as.character(cut2(bdE4$f.3710.0.0)))))

# Menstruating today #
levels(bdE4$f.3720.0.0) <- c(levels(bdE4$f.3720.0.0),"Menopause")
bdE4$f.3720.0.0[(bdE4$f.2724.0.0 != "Yes" & bdE4$f.2724.0.0 != "Not sure - had a hysterectomy")] <- "Menopause"

# Former alchol drinker #
levels(bdE4$f.3731.0.0) <- c(levels(bdE4$f.3731.0.0),"Drink alchol")
bdE4$f.3731.0.0[bdE4$f.1558.0.0 != "Never"] <- "Drink alchol"

# Stomach pain #
levels(bdE4$f.3741.0.0) <- c(levels(bdE4$f.3741.0.0),"No Stomach or abdominal pain")
bdE4$f.3741.0.0[bdE4$f.6159.0.l6 == 0] <- "No Stomach or abdominal pain"

# Chest pain when walking uphill #
levels(bdE4$f.3751.0.0) <- c(levels(bdE4$f.3751.0.0),"No chest pain")
bdE4$f.3751.0.0[bdE4$f.2335.0.0 != "Yes" | bdE4$f.3606.0.0=="Yes" ] <- "No chest pain"

# Age hay fever etc. #
bdE4$f.3761.0.0 <- ifelse(bdE4$f.3761.0.0==-1, NA, ifelse(bdE4$f.3761.0.0==-3,NA,bdE4$f.3761.0.0))
bdE4$f.3761.0.0 <- as.factor(ifelse(bdE4$f.6152.0.l6=="Hayfever, allergic rhinitis or eczema" ,as.character(cut2(bdE4$f.3761.0.0)),"No hay fever etc."))

# Knee pain #
levels(bdE4$f.3773.0.0) <- c(levels(bdE4$f.3773.0.0),"No knee pain")
bdE4$f.3773.0.0[bdE4$f.6159.0.l8 == 0] <- "No knee pain"

# Age Asthma #
bdE4$f.3786.0.0 <- ifelse(bdE4$f.3786.0.0==-1, NA, ifelse(bdE4$f.3786.0.0==-3,NA,bdE4$f.3786.0.0))
bdE4$f.3786.0.0 <- as.factor(ifelse(bdE4$f.6152.0.l5=="Asthma" ,as.character(cut2(bdE4$f.3786.0.0)),"No Asthma"))

# Headaches #
levels(bdE4$f.3799.0.0) <- c(levels(bdE4$f.3799.0.0),"No Headaches")
bdE4$f.3799.0.0[bdE4$f.6159.0.l2 == 0] <- "No Headaches"

# Time since PSA #
bdE4$f.3809.0.0 <- ifelse(bdE4$f.3809.0.0==-1, NA, ifelse(bdE4$f.3809.0.0==-3,NA,ifelse(bdE4$f.3809.0.0==-10,0,bdE4$f.3809.0.0)))
bdE4$f.3809.0.0 <- as.factor(ifelse(bdE4$f.2365.0.0=="Yes" ,as.character(cut2(bdE4$f.3809.0.0)),"No PSA"))

# Number of stillbirths - make two categories 0, >0 #
bdE4$f.3829.0.0 <- as.factor(ifelse(bdE4$f.3829.0.0==-1, NA, ifelse(bdE4$f.3829.0.0==-3,NA,ifelse(bdE4$f.3829.0.0==0,1,2))))
levels(bdE4$f.3829.0.0) <- c(levels(bdE4$f.3829.0.0),"Not stillbirth")
bdE4$f.3829.0.0[bdE4$f.2774.0.0!="Yes"] <- "Not stillbirth"

# Number of spontaneous miscarriage #
bdE4$f.3839.0.0 <- ifelse(bdE4$f.3839.0.0==-1, NA, ifelse(bdE4$f.3839.0.0==-3,NA,bdE4$f.3839.0.0))
bdE4$f.3839.0.0 <- as.factor(ifelse(bdE4$f.2774.0.0=="Yes" ,as.character(cut2(bdE4$f.3839.0.0)),"No miscarriage"))

# Number of pregnancy termination #
bdE4$f.3849.0.0 <- ifelse(bdE4$f.3849.0.0==-1, NA, ifelse(bdE4$f.3849.0.0==-3,NA,bdE4$f.3849.0.0))
bdE4$f.3849.0.0 <- as.factor(ifelse(bdE4$f.2774.0.0=="Yes" ,as.character(cut2(bdE4$f.3849.0.0)),"No miscarriage"))

# Reason to stop drinking alchol #
levels(bdE4$f.3859.0.0) <- c(levels(bdE4$f.3859.0.0),"Drink alchol or no previous alchol")
bdE4$f.3859.0.0[bdE4$f.1558.0.0 != "Never" | bdE4$f.3731.0.0!="Yes"] <- "Drink alchol or no previous alchol"

# Age of primiparous women at birth of child #
bdE4$f.3872.0.0 <- ifelse(bdE4$f.3872.0.0==-3,NA,ifelse(bdE4$f.3872.0.0==-4,NA,bdE4$f.3872.0.0))
bdE4$f.3872.0.0 <- as.factor(ifelse(bdE4$f.2734.0.0b == 1  ,as.character(cut2(bdE4$f.3872.0.0)),"No primiparous"))

# Age at bilateral oophorectomy #
bdE4$f.3882.0.0 <- ifelse(bdE4$f.3882.0.0==-1, NA, ifelse(bdE4$f.3882.0.0==-3,NA,bdE4$f.3882.0.0))
bdE4$f.3882.0.0 <- as.factor(ifelse(bdE4$f.2834.0.0 == "Yes" ,as.character(cut2(bdE4$f.3882.0.0)),"No primiparous"))

# Age heart attack #
bdE4$f.3894.0.0 <- ifelse(bdE4$f.3894.0.0==-1, NA, ifelse(bdE4$f.3894.0.0==-3,NA,bdE4$f.3894.0.0))
bdE4$f.3894.0.0 <- as.factor(ifelse(bdE4$f.6150.0.l2 == "Heart attack" ,as.character(cut2(bdE4$f.3894.0.0)),"No heart attack"))

# Adopted father still alive #
levels(bdE4$f.3912.0.0) <- c(levels(bdE4$f.3912.0.0),"No adopted")
bdE4$f.3912.0.0[bdE4$f.1767.0.0 != "Yes"] <- "No adopted"

# Adopted mother still alive #
levels(bdE4$f.3942.0.0) <- c(levels(bdE4$f.3942.0.0),"No adopted")
bdE4$f.3942.0.0[bdE4$f.1767.0.0 != "Yes"] <- "No adopted"

# N. adopted brother #
bdE4$f.3972.0.0b <- bdE4$f.3972.0.0
bdE4$f.3972.0.0 <- ifelse(bdE4$f.3972.0.0==-1, NA, ifelse(bdE4$f.3972.0.0==-3,NA,bdE4$f.3972.0.0))
bdE4$f.3972.0.0 <- as.factor(ifelse(bdE4$f.1767.0.0 == "Yes" ,as.character(cut2(bdE4$f.3972.0.0)),"No adopted"))

# N. adopted sisters #
bdE4$f.3982.0.0b <- bdE4$f.3982.0.0
bdE4$f.3982.0.0 <- ifelse(bdE4$f.3982.0.0==-1, NA, ifelse(bdE4$f.3982.0.0==-3,NA,bdE4$f.3982.0.0))
bdE4$f.3982.0.0 <- as.factor(ifelse(bdE4$f.1767.0.0 == "Yes" ,as.character(cut2(bdE4$f.3982.0.0)),"No adopted"))

# Age emphysema/chronic bronchitis #
bdE4$f.3992.0.0 <- ifelse(bdE4$f.3992.0.0==-1, NA, ifelse(bdE4$f.3992.0.0==-3,NA,bdE4$f.3992.0.0))
bdE4$f.3992.0.0 <- as.factor(ifelse(bdE4$f.6152.0.l3=="Emphysema/chronic bronchitis" ,as.character(cut2(bdE4$f.3992.0.0)),"No emphysema etc."))

# Age deep-vein thrombosis #
bdE4$f.4012.0.0 <- ifelse(bdE4$f.4012.0.0==-1, NA, ifelse(bdE4$f.4012.0.0==-3,NA,bdE4$f.4012.0.0))
bdE4$f.4012.0.0 <- as.factor(ifelse(bdE4$f.6152.0.l2=="Blood clot in the leg (DVT)" ,as.character(cut2(bdE4$f.4012.0.0)),"No deep-vein thrombosis"))

# Age pulmunary embolism #
bdE4$f.4022.0.0 <- ifelse(bdE4$f.4022.0.0==-1, NA, ifelse(bdE4$f.4022.0.0==-3,NA,bdE4$f.4022.0.0))
bdE4$f.4022.0.0 <- as.factor(ifelse(bdE4$f.6152.0.l4=="Blood clot in the lung" ,as.character(cut2(bdE4$f.4022.0.0)),"No pulmunary embolism"))

# Gestational diabetes #
levels(bdE4$f.4041.0.0) <- c(levels(bdE4$f.4041.0.0),"No diabetes")
bdE4$f.4041.0.0[bdE4$f.2443.0.0 != "Yes"] <- "No diabetes"

# Age stroke diagnosed #
bdE4$f.4056.0.0 <- ifelse(bdE4$f.4056.0.0==-1, NA, ifelse(bdE4$f.4056.0.0==-3,NA,bdE4$f.4056.0.0))
bdE4$f.4056.0.0 <- as.factor(ifelse(bdE4$f.6150.0.l4 == "Stroke" ,as.character(cut2(bdE4$f.4056.0.0)),"No Stroke"))

# Facial pain #
levels(bdE4$f.4067.0.0) <- c(levels(bdE4$f.4067.0.0),"No Facial pain")
bdE4$f.4067.0.0[bdE4$f.6159.0.l3 == 0] <- "No Facial pain"

# Diastolic  and systolic blood pressure - use manual reading if missing and average #
bdE4$f.4079.0.0 <- ifelse(is.na(bdE4$f.4079.0.0),bdE4$f.94.0.0,bdE4$f.4079.0.0)
bdE4$f.4079.0.1 <- ifelse(is.na(bdE4$f.4079.0.1),bdE4$f.94.0.1,bdE4$f.4079.0.1)
bdE4$f.4079.0.0  <- rowMeans(cbind(bdE4$f.4079.0.0,bdE4$f.4079.0.1), na.rm=T)
bdE4$f.4079.0.1 <- NULL
bdE4$f.4080.0.0 <- ifelse(is.na(bdE4$f.4080.0.0),bdE4$f.93.0.0,bdE4$f.4080.0.0)
bdE4$f.4080.0.1 <- ifelse(is.na(bdE4$f.4080.0.1),bdE4$f.93.0.1,bdE4$f.4080.0.1)
bdE4$f.4080.0.0  <- rowMeans(cbind(bdE4$f.4080.0.0,bdE4$f.4080.0.1), na.rm=T)
bdE4$f.4080.0.1 <- NULL
bdE4$f.93.0.0 <- NULL
bdE4$f.93.0.1 <- NULL
bdE4$f.94.0.0 <- NULL
bdE4$f.94.0.1 <- NULL

print("four")

# Method assignment blood pressure #
bdE4$f.4081.0.0 <- NULL
bdE4$f.4081.0.1 <- NULL

# Heel ultrasound method #
bdE4$f.4092.0.0 <- NULL
bdE4$f.4095.0.0 <- NULL

# Fractured heel #
bdE4$f.4093.0.0 <- NULL
bdE4$f.4096.0.0 <- NULL

# Exclude these manual masured bone densities because not clear if left or right #
bdE4$f.78.0.0 <- NULL
bdE4$f.77.0.0 <- NULL

# Heel bone density - set to missing because left and right specific (and we already have the overall measurement) #
bdE4$f.4100.0.0 <- NULL
bdE4$f.4101.0.0 <- NULL
bdE4$f.4104.0.0 <- NULL
bdE4$f.4106.0.0 <- NULL
bdE4$f.4119.0.0 <- NULL
bdE4$f.4120.0.0 <- NULL
bdE4$f.4123.0.0 <- NULL
bdE4$f.4105.0.0 <- NULL
bdE4$f.4124.0.0 <- NULL
bdE4$f.4106.0.0 <- NULL
bdE4$f.4125.0.0 <- NULL
bdE4$f.4141.0.0 <- NULL
bdE4$f.4146.0.0 <- NULL
bdE4$f.4139.0.0 <- NULL
bdE4$f.4144.0.0 <- NULL
bdE4$f.4140.0.0 <- NULL
bdE4$f.4145.0.0 <- NULL
bdE4$f.4138.0.0 <- NULL
bdE4$f.4143.0.0 <- NULL

# Arterial pulse-wave stiffness device ID #
bdE4$f.4136.0.0 <- NULL

# Stifness method #
bdE4$f.4186.0.0 <- NULL

# Arterial stiffness device ID #
bdE4$f.4206.0.0 <- NULL

# Pulse wave velocity (manual) #
bdE4$f.4207.0.0 <- NULL

# EXCLUDE ALL THE HEARING TEST VARIABLES WHICH ARE NOT USEFUL #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("4229","4230","4233","4234","4235","4236","4237","4238","4239","4240","4241","4242","4243","4244","4245","4246","4247","4248","4249","4252","4269","4276","4270","4277","4272","4279")]

# Total triplet correct #
bdE4$f.4232.0.M <- rowSums(cbind(bdE4$f.4232.0.1=="yes",
bdE4$f.4232.0.2=="yes",
bdE4$f.4232.0.3=="yes",
bdE4$f.4232.0.4=="yes",
bdE4$f.4232.0.5=="yes",
bdE4$f.4232.0.6=="yes",
bdE4$f.4232.0.7=="yes",
bdE4$f.4232.0.8=="yes",
bdE4$f.4232.0.9=="yes",
bdE4$f.4232.0.10=="yes",
bdE4$f.4232.0.11=="yes",
bdE4$f.4232.0.12=="yes",
bdE4$f.4232.0.13=="yes",
bdE4$f.4232.0.14=="yes",
bdE4$f.4232.0.15=="yes"), na.rm=F)
bdE4$f.4232.0.1 <- NULL
bdE4$f.4232.0.2 <- NULL
bdE4$f.4232.0.3 <- NULL
bdE4$f.4232.0.4 <- NULL
bdE4$f.4232.0.5 <- NULL
bdE4$f.4232.0.6 <- NULL
bdE4$f.4232.0.7 <- NULL
bdE4$f.4232.0.8 <- NULL
bdE4$f.4232.0.9 <- NULL
bdE4$f.4232.0.10 <- NULL
bdE4$f.4232.0.11 <- NULL
bdE4$f.4232.0.12 <- NULL
bdE4$f.4232.0.13 <- NULL
bdE4$f.4232.0.14 <- NULL
bdE4$f.4232.0.15 <- NULL

# EXCLUDE THE MEMORY TEST VARIABLES WHICH ARE NOT USEFULL #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("4250","4251","4253","4254","4255","4256","4257","4258","4259","4260","4268")]

# Completion status #
bdE4$f.4275.0.0 <- NULL
bdE4$f.4281.0.0 <- NULL
bdE4$f.4287.0.0 <- NULL

# maximum digits rememebered #
bdE4$f.4282.0.0 <- ifelse(bdE4$f.4282.0.0==-1, 0,bdE4$f.4282.0.0)

# Total time to complete test #
bdE4$f.4285.0.0 <- NULL

# Time initial screen shown #
bdE4$f.4286.0.0 <- NULL 

# Time to answer #
bdE4$f.4288.0.0 <- NULL 

# Time screen exit #
bdE4$f.4289.0.0 <- NULL 

# Duration screen displayed #
bdE4$f.4290.0.0 <- NULL 

# Initial answer #
bdE4$f.4292.0.0 <- NULL

# Final answer #
bdE4$f.4293.0.0 <- NULL

# History of attempts #
bdE4$f.4295.0.0 <- NULL

# Monthly wine take #
bdE4$f.4407.0.0 <- ifelse(bdE4$f.4407.0.0==-1, NA, ifelse(bdE4$f.4407.0.0==-3,NA,bdE4$f.4407.0.0))
bdE4$f.4407.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Special occasions only" | bdE4$f.1558.0.0=="One to three times a month",as.character(cut2(bdE4$f.4407.0.0)),"Not drinking rarerly"))

# Monthly white wine #
bdE4$f.4418.0.0 <- ifelse(bdE4$f.4418.0.0==-1, NA, ifelse(bdE4$f.4418.0.0==-3,NA,bdE4$f.4418.0.0))
bdE4$f.4418.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Special occasions only" | bdE4$f.1558.0.0=="One to three times a month",as.character(cut2(bdE4$f.4418.0.0)),"Not drinking rarerly"))

# Monthly beer + cider #
bdE4$f.4429.0.0 <- ifelse(bdE4$f.4429.0.0==-1, NA, ifelse(bdE4$f.4429.0.0==-3,NA,bdE4$f.4429.0.0))
bdE4$f.4429.0.0 <- as.factor(ifelse(bdE4$f.1558.0.0=="Special occasions only" | bdE4$f.1558.0.0=="One to three times a month",as.character(cut2(bdE4$f.4429.0.0)),"Not drinking rarerly"))

# Monthly spirits - make two categories 0, >0  #
bdE4$f.4440.0.0 <-  as.factor(ifelse(bdE4$f.4440.0.0==-1, NA, ifelse(bdE4$f.4440.0.0==-3,NA,ifelse(bdE4$f.4440.0.0==0,0,">=1"))))
levels(bdE4$f.4440.0.0) <- c(levels(bdE4$f.4440.0.0),"Not drinking rarerly")
bdE4$f.4440.0.0[bdE4$f.1558.0.0!="Special occasions only" & bdE4$f.1558.0.0!="One to three times a month"] <- "Not drinking rarerly"

# Monthly fortified wine - make two categories 0, >0 #
bdE4$f.4451.0.0 <-  as.factor(ifelse(bdE4$f.4451.0.0==-1, NA, ifelse(bdE4$f.4451.0.0==-3,NA,ifelse(bdE4$f.4451.0.0==0,0,">=1"))))
levels(bdE4$f.4451.0.0) <- c(levels(bdE4$f.4451.0.0),"Not drinking rarerly")
bdE4$f.4451.0.0[bdE4$f.1558.0.0!="Special occasions only" & bdE4$f.1558.0.0!="One to three times a month"] <- "Not drinking rarerly"

# Monthly other wine - make two categories 0, >0 #
bdE4$f.4462.0.0 <-  as.factor(ifelse(bdE4$f.4462.0.0==-1, NA, ifelse(bdE4$f.4462.0.0==-3,NA,ifelse(bdE4$f.4462.0.0==0,0,">=1"))))
levels(bdE4$f.4462.0.0) <- c(levels(bdE4$f.4462.0.0),"Not drinking rarerly")
bdE4$f.4462.0.0[bdE4$f.1558.0.0!="Special occasions only" & bdE4$f.1558.0.0!="One to three times a month"] <- "Not drinking rarerly"

# Longest depression period #
bdE4$f.4609.0.0 <- ifelse(bdE4$f.4609.0.0==-1, NA, ifelse(bdE4$f.4609.0.0==-3,NA,bdE4$f.4609.0.0))
bdE4$f.4609.0.0 <- as.factor(ifelse(bdE4$f.4598.0.0=="Yes",as.character(cut2(bdE4$f.4609.0.0)),"No depressed"))

# N. depression episodes #
bdE4$f.4620.0.0 <- ifelse(bdE4$f.4620.0.0==-1, NA, ifelse(bdE4$f.4620.0.0==-3,NA,bdE4$f.4620.0.0))
bdE4$f.4620.0.0 <- as.factor(ifelse(bdE4$f.4598.0.0=="Yes",as.character(cut2(bdE4$f.4620.0.0)),"No depressed"))

# Age glaucoma #
bdE4$f.4689.0.0 <- ifelse(bdE4$f.4689.0.0==-1, NA, ifelse(bdE4$f.4689.0.0==-3,NA,bdE4$f.4689.0.0))
bdE4$f.4689.0.0 <- as.factor(ifelse(bdE4$f.6148.0.l3=="Glaucoma" ,as.character(cut2(bdE4$f.4689.0.0)),"No Glaucoma"))

# Age cataract #
bdE4$f.4700.0.0 <- ifelse(bdE4$f.4700.0.0==-1, NA, ifelse(bdE4$f.4700.0.0==-3,NA,bdE4$f.4700.0.0))
bdE4$f.4700.0.0 <- as.factor(ifelse(bdE4$f.6148.0.l5=="Cataract" ,as.character(cut2(bdE4$f.4700.0.0)),"No Cataract"))

# Tinnitus severity #
levels(bdE4$f.4814.0.0) <- c(levels(bdE4$f.4814.0.0),"No tinnitus")
bdE4$f.4814.0.0[bdE4$f.4803.0.0!="Yes, now most or all of the time" & bdE4$f.4803.0.0!="Yes, now a lot of the time" & bdE4$f.4803.0.0!="Yes, now some of the time" & bdE4$f.4803.0.0!="Yes, but not now, but have in the past"] <- "No tinnitus"


# FLUID INTELLIGENT TEST VARIABLES TO EXCLUDE #
bdE4$f.4935.0.0 <- NULL
bdE4$f.4946.0.0 <- NULL
bdE4$f.4957.0.0 <- NULL
bdE4$f.4968.0.0 <- NULL
bdE4$f.4979.0.0 <- NULL
bdE4$f.4990.0.0 <- NULL
bdE4$f.5001.0.0 <- NULL
bdE4$f.5012.0.0 <- NULL
bdE4$f.5556.0.0 <- NULL
bdE4$f.5699.0.0 <- NULL
bdE4$f.5779.0.0 <- NULL
bdE4$f.5790.0.0 <- NULL
bdE4$f.5866.0.0 <- NULL

# N. of older sibilings #
bdE4$f.5057.0.0 <- ifelse(bdE4$f.5057.0.0==-1, NA, ifelse(bdE4$f.5057.0.0==-3,NA,bdE4$f.5057.0.0))
bdE4$f.5057.0.0 <- as.factor(ifelse(bdE4$f.1767.0.0=="Yes","Adopted",as.character(cut2(bdE4$f.5057.0.0))))

# VISUAL ACUITY TEST VARIABLES TO EXCLUDE #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("5074","5075","5078","5079","5080","5081","5082","5083","5186","5187","5188","5189","5189","5191","5194","5196","5199","5200","5201","5204","5205","5206","5207","5208","5211","5212")]

# Keep n. correct in round 1, exclude the other #
bdE4$f.5076.0.M <- rowSums(cbind(bdE4$f.5076.0.0,
bdE4$f.5076.0.1,
bdE4$f.5076.0.2,
bdE4$f.5076.0.3,
bdE4$f.5076.0.4,
bdE4$f.5076.0.5,
bdE4$f.5076.0.6,
bdE4$f.5076.0.7,
bdE4$f.5076.0.8,
bdE4$f.5076.0.9,
bdE4$f.5076.0.10,
bdE4$f.5076.0.11,
bdE4$f.5076.0.12,
bdE4$f.5076.0.13,
bdE4$f.5076.0.14,
bdE4$f.5076.0.15), na.rm=T)
bdE4$f.5076.0.M[bdE4$f.5076.0.M==0] <- NA

bdE4$f.5077.0.M <- rowSums(cbind(bdE4$f.5077.0.0,
bdE4$f.5077.0.1,
bdE4$f.5077.0.2,
bdE4$f.5077.0.3,
bdE4$f.5077.0.4,
bdE4$f.5077.0.5,
bdE4$f.5077.0.6,
bdE4$f.5077.0.7,
bdE4$f.5077.0.8,
bdE4$f.5077.0.9,
bdE4$f.5077.0.10,
bdE4$f.5077.0.11,
bdE4$f.5077.0.12,
bdE4$f.5077.0.13,
bdE4$f.5077.0.14,
bdE4$f.5077.0.15), na.rm=T)
bdE4$f.5077.0.M[bdE4$f.5077.0.M==0] <- NA

bdE4$f.5077.0.0 <- NULL
bdE4$f.5077.0.1 <- NULL
bdE4$f.5077.0.2 <- NULL
bdE4$f.5077.0.3 <- NULL
bdE4$f.5077.0.4 <- NULL
bdE4$f.5077.0.5 <- NULL
bdE4$f.5077.0.6 <- NULL
bdE4$f.5077.0.7 <- NULL
bdE4$f.5077.0.8 <- NULL
bdE4$f.5077.0.9 <- NULL
bdE4$f.5077.0.10 <- NULL
bdE4$f.5077.0.11 <- NULL
bdE4$f.5077.0.12 <- NULL
bdE4$f.5077.0.13 <- NULL
bdE4$f.5077.0.14 <- NULL
bdE4$f.5077.0.15 <- NULL
bdE4$f.5076.0.0 <- NULL
bdE4$f.5076.0.1 <- NULL
bdE4$f.5076.0.2 <- NULL
bdE4$f.5076.0.3 <- NULL
bdE4$f.5076.0.4 <- NULL
bdE4$f.5076.0.5 <- NULL
bdE4$f.5076.0.6 <- NULL
bdE4$f.5076.0.7 <- NULL
bdE4$f.5076.0.8 <- NULL
bdE4$f.5076.0.9 <- NULL
bdE4$f.5076.0.10 <- NULL
bdE4$f.5076.0.11 <- NULL
bdE4$f.5076.0.12 <- NULL
bdE4$f.5076.0.13 <- NULL
bdE4$f.5076.0.14 <- NULL
bdE4$f.5076.0.15 <- NULL


# REFRACTOMETRY RESULTS, KEEP THE SECOND ROUND BECAUSE OVERALL THE BEST (see f.5221.0.0) #
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("5084","5085","5086","5087","5088","5089","5096","5097","5098","5099","5100","5101","5102","5103","5104","5105","5106","5107","5108","5109","5110","5111","5112","5113","5114","5115","5116","5117","5118","5119","5132","5133","5134","5135","5156","5157","5158","5159","5160","5161","5162","5163") & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4))%in%c("0","2","3","4","5","6","7","8","9"))]

# REFRACTOMETRY RESULTS, EXCLUDE ALL THESE VARIABLES ##
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("5090","5091","5136","5137","5138","5139","5140","5141","5142","5143","5144","5145","5146","5147","5148","5149","5152","5155","5164","5185","5198","5215","5221","5237","5251","5274","5276","5292","5306","5324","5325","5326","5327","5328")]

# Average week intake of other alchool - make two categories 0, >0 #
bdE4$f.5364.0.0 <- as.factor(ifelse(bdE4$f.5364.0.0==-1, NA, ifelse(bdE4$f.5364.0.0==-3,NA,ifelse(bdE4$f.5364.0.0==0,0,">= 1"))))
levels(bdE4$f.5364.0.0) <- c(levels(bdE4$f.5364.0.0),"Not drinking / Not drinking often")
bdE4$f.5364.0.0[bdE4$f.1558.0.0!="Daily or almost daily" & bdE4$f.1558.0.0!="Three or four times a week" & bdE4$f.1558.0.0!="Once or twice a week"] <- "Not drinking / Not drinking often"

# Longest period of unenthusiasm #
bdE4$f.5375.0.0 <- ifelse(bdE4$f.5375.0.0==-1, NA, ifelse(bdE4$f.5375.0.0==-3,NA,bdE4$f.5375.0.0))
bdE4$f.5375.0.0 <- as.factor(ifelse(bdE4$f.4631.0.0=="Yes" ,as.character(cut2(bdE4$f.5375.0.0)),"Yes Enthusiam"))

# N. unenthusiasm periods #
bdE4$f.5386.0.0 <- ifelse(bdE4$f.5386.0.0==-1, NA, ifelse(bdE4$f.5386.0.0==-3,NA,bdE4$f.5386.0.0))
bdE4$f.5386.0.0 <- as.factor(ifelse(bdE4$f.4631.0.0=="Yes" ,as.character(cut2(bdE4$f.5386.0.0)),"Yes Enthusiam"))

# Eye affected by Ambloyopia #
levels(bdE4$f.5408.0.0) <- c(levels(bdE4$f.5408.0.0),"No lazy eye")
bdE4$f.5408.0.0[bdE4$f.6147.0.l6 == 0] <- "No lazy eye"

# Eye affected by trauma #
levels(bdE4$f.5419.0.0) <- c(levels(bdE4$f.5419.0.0),"No trauma")
bdE4$f.5419.0.0[bdE4$f.6148.0.l4 == 0] <- "No trauma"

# Age lost of vision - deleted because too less subjects (< 20 ) #
bdE4$f.5430.0.0 <- NULL

# Eye affected by cataract #
levels(bdE4$f.5441.0.0) <- c(levels(bdE4$f.5441.0.0),"No cataract")
bdE4$f.5441.0.0[bdE4$f.6148.0.l5 == 0 ] <- "No cataract"

# Leg pain when standing still #
levels(bdE4$f.5452.0.0) <- c(levels(bdE4$f.5452.0.0),"No leg pain")
bdE4$f.5452.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Leg pain in claf #
levels(bdE4$f.5463.0.0) <- c(levels(bdE4$f.5463.0.0),"No leg pain")
bdE4$f.5463.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Leg pain uphill #
levels(bdE4$f.5474.0.0) <- c(levels(bdE4$f.5474.0.0),"No leg pain")
bdE4$f.5474.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Leg normally walking #
levels(bdE4$f.5485.0.0) <- c(levels(bdE4$f.5485.0.0),"No leg pain")
bdE4$f.5485.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Leg pain disappears #
levels(bdE4$f.5496.0.0) <- c(levels(bdE4$f.5496.0.0),"No leg pain")
bdE4$f.5496.0.0[bdE4$f.4728.0.0 != "Yes" | bdE4$f.5485.0.0 !="Yes"] <- "No leg pain"

# Leg pain action taken #
levels(bdE4$f.5507.0.0) <- c(levels(bdE4$f.5507.0.0),"No leg pain")
bdE4$f.5507.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Leg pain effect on standing still #
levels(bdE4$f.5518.0.0) <- c(levels(bdE4$f.5518.0.0),"No leg pain")
bdE4$f.5518.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Surgery on leg arteries #
levels(bdE4$f.5529.0.0) <- c(levels(bdE4$f.5529.0.0),"No leg pain")
bdE4$f.5529.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Surgery/ amputation #
levels(bdE4$f.5540.0.0) <- c(levels(bdE4$f.5540.0.0),"No leg pain")
bdE4$f.5540.0.0[bdE4$f.4728.0.0 != "Yes" ] <- "No leg pain"

# Which eye presbopia #
levels(bdE4$f.5610.0.0) <- c(levels(bdE4$f.5610.0.0),"No presbyopia")
bdE4$f.5610.0.0[bdE4$f.6147.0.l3 == 0] <- "No presbyopia"

# Longest manic episode #
levels(bdE4$f.5663.0.0) <- c(levels(bdE4$f.5663.0.0),"I'm not manic/irritable")
bdE4$f.5663.0.0[bdE4$f.4642.0.0 != "Yes"  & bdE4$f.4653.0.0!="Yes"] <- "I'm not manic/irritable"

print("five")

# Severity of manic #
levels(bdE4$f.5674.0.0) <- c(levels(bdE4$f.5674.0.0),"I'm not manic/irritable")
bdE4$f.5674.0.0[bdE4$f.4642.0.0 != "Yes" & bdE4$f.4653.0.0!="Yes"] <- "I'm not manic/irritable"

# Which eye hypermetropia #
levels(bdE4$f.5832.0.0) <- c(levels(bdE4$f.5832.0.0),"No hypermetropia")
bdE4$f.5832.0.0[bdE4$f.6147.0.l2 == 0] <- "No hypermetropia"

# Which eye myopia #
levels(bdE4$f.5843.0.0) <- c(levels(bdE4$f.5843.0.0),"No myopia")
bdE4$f.5843.0.0[bdE4$f.6147.0.l1 == 0] <- "No myopia"

# Which eye astigmatism #
levels(bdE4$f.5855.0.0) <- c(levels(bdE4$f.5855.0.0),"No astigmatism")
bdE4$f.5855.0.0[bdE4$f.6147.0.l4 == 0 ] <- "No astigmatism"

# Which eye other conditions #
levels(bdE4$f.5877.0.0) <- c(levels(bdE4$f.5877.0.0),"No other")
bdE4$f.5877.0.0[bdE4$f.6147.0.l7 == 0 ] <- "No other"

# Which eye diabetes #
levels(bdE4$f.5890.0.0) <- c(levels(bdE4$f.5890.0.0),"No diab eye")
bdE4$f.5890.0.0[bdE4$f.6148.0.l2 == 0] <- "No diab eye"

# Age eye diabetes #
bdE4$f.5901.0.0 <- ifelse(bdE4$f.5901.0.0==-1, NA, ifelse(bdE4$f.5901.0.0==-3,NA,bdE4$f.5901.0.0))
bdE4$f.5901.0.0 <- as.factor(ifelse(bdE4$f.6148.0.l2 == "Diabetes related eye disease" ,as.character(cut2(bdE4$f.5901.0.0)),"No diab eye"))

# Which eye macular #
levels(bdE4$f.5912.0.0) <- c(levels(bdE4$f.5912.0.0),"No macular")
bdE4$f.5912.0.0[bdE4$f.6148.0.l6 == 0 ] <- "No macular"

# Age macular - too less events #
bdE4$f.5923.0.0 <- NULL

# Which eye other #
levels(bdE4$f.5934.0.0) <- c(levels(bdE4$f.5934.0.0),"No eye other")
bdE4$f.5934.0.0[bdE4$f.6148.0.l7 == 0 ] <- "No eye other"

# Age eye other # 
bdE4$f.5945.0.0 <- ifelse(bdE4$f.5945.0.0==-1, NA, ifelse(bdE4$f.5945.0.0==-3,NA,bdE4$f.5945.0.0)) 
bdE4$f.5945.0.0 <- as.factor(ifelse(bdE4$f.6148.0.l7 == "Other serious eye condition",as.character(cut2(bdE4$f.5945.0.0)),"No eye other"))

# Previsouly smoked cigarettes all days #
levels(bdE4$f.5959.0.0) <- c(levels(bdE4$f.5959.0.0),"No smoke pipe etc")
bdE4$f.5959.0.0[bdE4$f.1239.0.0 != "Yes, on most or all days" | bdE4$f.3446.0.0!="Cigars or pipes" ] <- "No smoke pipe etc"

# ECG VARIABLES TO EXCLUDE, TOO MANY MISSING (>80%) #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("5984","5985","5986","5988","5991","5992","5993","6023","6038","6039","5983","5987","6014","6015","6016","6017","6019","6020","6022","6024","6032","6033","6034")]

# OCT measures #
bdE4$f.6070.0.0 <- NULL
bdE4$f.6072.0.0 <- NULL

# Glass required #
bdE4$f.6074.0.0 <- NULL
bdE4$f.6075.0.0 <- NULL

# Which eye glaucoma #
levels(bdE4$f.6119.0.0) <- c(levels(bdE4$f.6119.0.0),"No Glaucoma")
bdE4$f.6119.0.0[bdE4$f.6148.0.l3 == 0 ] <- "No Glaucoma"

# N. of cigarettes previsouly smoked - remove because no subjects #
bdE4$f.6183.0.0 <- NULL

# Age stop smoking - remove because no subjects #
bdE4$f.6194.0.0 <- NULL

# Which eye(s) affected by strabismus (squint) #
levels(bdE4$f.6205.0.0) <- c(levels(bdE4$f.6205.0.0),"No strabismus")
bdE4$f.6205.0.0[bdE4$f.6147.0.l5 == 0] <- "No strabismus"

# THESE VARIABLES EXCLUDED BECAUSE ALREADY SUMMARIZED BY 134, 135, 136, 137 #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("20001","20002","20003","20004","20006")]

# Age cancer diagnosed - Problems in the analysis, excluded #
bdE4$f.20007.0.0 <- NULL

# Age non-cancer illness #
bdE4$f.20009.0.0 <- ifelse(bdE4$f.20009.0.0==-1, NA, ifelse(bdE4$f.20009.0.0==-2,NA,ifelse(bdE4$f.20009.0.0==-3,NA,bdE4$f.20009.0.0)))
bdE4$f.20009.0.0 <- as.factor(ifelse(bdE4$f.135.0.0b > 0  ,as.character(cut2(bdE4$f.20009.0.0)),"No illness"))

# Age operation #
bdE4$f.20011.0.0 <- ifelse(bdE4$f.20011.0.0==-1, NA, ifelse(bdE4$f.20011.0.0==-2,NA,ifelse(bdE4$f.20011.0.0==-3,NA,bdE4$f.20011.0.0)))
bdE4$f.20011.0.0 <- as.factor(ifelse(bdE4$f.136.0.0b > 0  ,as.character(cut2(bdE4$f.20011.0.0)),"No operation"))

# AGE AT MEDICAL CONDITIONS - REMOVED ALL THE OTHER ARRAYS EXCEPT FIRST #
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("20007","20009","20011") & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(1,31,1)))]

# MEDICAL CONDITIONS - REMOVED ALL THESE VARIABLES #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("20008","20010","20012","20013","20014")]

# Home area population density - urban or rural #
# Recode because too many levels #
bdE4$f.20118.0.0 <- as.factor(ifelse(bdE4$f.20118.0.0%in%c("England/Wales - Urban - sparse","England/Wales - Urban - less sparse","Scotland - Large Urban Area","Scotland - Other Urban Area"),"Urban",
ifelse(bdE4$f.20118.0.0%in%c("England/Wales - Town and Fringe - sparse","England/Wales - Town and Fringe - less sparse","Scotland - Accessible Small Town","Scotland - Remote Small Town","Scotland - Very Remote Small Town"),"Town",
ifelse(bdE4$f.20118.0.0%in%c("England/Wales - Village - sparse","England/Wales - Village - less sparse"),"Village",
ifelse(bdE4$f.20118.0.0%in%c("England/Wales - Hamlet and Isolated dwelling - sparse","England/Wales - Hamlet and Isolated Dwelling - less sparse","Scotland - Accessible Rural","Scotland - Remote Rural","Scotland - Very Remote Rural"),"Rural",NA)))))

# Ethnic background #
# Recode because too many levels #
bdE4$f.21000.0.0 <- as.factor(ifelse(bdE4$f.21000.0.0%in%c("White","Any other white background","British","Irish"),"White",
ifelse(bdE4$f.21000.0.0%in%c("Black or Black British","Any other Black background","African","Caribbean"),"Black",
ifelse(bdE4$f.21000.0.0%in%c("Indian","Pakistani","Bangladeshi","Any other Asian background","Asian or Asian British"),"Asian",
ifelse(bdE4$f.21000.0.0%in%c("Chinese"), "Chinese",
ifelse(bdE4$f.21000.0.0%in%c("Mixed","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background"), "Mixed",
ifelse(bdE4$f.21000.0.0%in%c("Other ethnic group"), "Other", NA)))))))


# Acceptability of each blow result (text) #
bdE4$f.20031.0.0 <- NULL
bdE4$f.20031.0.1 <- NULL
bdE4$f.20031.0.2 <- NULL

# Reasons for skipping #
bdE4$f.20041.0.0 <- NULL
bdE4$f.20042.0.0 <- NULL
bdE4$f.20043.0.0 <- NULL
bdE4$f.20044.0.0 <- NULL
bdE4$f.20045.0.0 <- NULL
bdE4$f.20046.0.0 <- NULL
bdE4$f.20047.0.0 <- NULL
bdE4$f.20048.0.0 <- NULL
bdE4$f.20051.0.0 <- NULL
bdE4$f.20052.0.0 <- NULL
bdE4$f.20053.0.0 <- NULL
bdE4$f.20054.0.0 <- NULL
bdE4$f.20055.0.0 <- NULL
bdE4$f.20056.0.0 <- NULL
bdE4$f.20057.0.0 <- NULL
bdE4$f.20058.0.0 <- NULL
bdE4$f.20059.0.0 <- NULL
bdE4$f.20060.0.0 <- NULL
bdE4$f.20061.0.0 <- NULL
bdE4$f.20062.0.0 <- NULL

# Illness of father #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="20107" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"Adopted")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.1767.0.0 == "Yes"] <- "Adopted"
		bdE4[paste("f.20107.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.20107.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="20107" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,10,1)))]
# Exclude these variables because "Don't know", "Prefer not to answer", etc.. #
bdE4$f.20107.0.l1 <- NULL
bdE4$f.20107.0.l2 <- NULL
bdE4$f.20107.0.l3 <- NULL
bdE4$f.20107.0.l4 <- NULL
bdE4$f.20107.0.l5 <- NULL
bdE4$f.20107.0.l6 <- NULL
bdE4$f.20107.0.l11 <- NULL

# Illness of mother #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="20110" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
		ll <- levels(tt[,1])[k]
		te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
		tef <- as.factor(te*1)
		levels(tef) <- c(levels(tef),"Adopted")
		tef[rowSums(is.na(tt))==ncol(tt)] <- NA
		tef[bdE4$f.1767.0.0 == "Yes"] <- "Adopted"
		bdE4[paste("f.20110.0.l",k,sep="")] <- tef
		levels(bdE4[[paste("f.20110.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="20110" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,12,1)))]
# Exclude these variables because "Don't know", "Prefer not to answer", etc.. #
bdE4$f.20110.0.l1 <- NULL
bdE4$f.20110.0.l2 <- NULL
bdE4$f.20110.0.l3 <- NULL
bdE4$f.20110.0.l4 <- NULL
bdE4$f.20110.0.l5 <- NULL
bdE4$f.20110.0.l6 <- NULL
bdE4$f.20110.0.l18 <- NULL

# Illness of sibilings #
tt <- bdE4[,unlist(lapply(strsplit(names(bdE4),"\\."),"[",2))=="20111" & !is.na(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)))]
for(k in 1:length(levels(tt[,1])))
	{
			ll <- levels(tt[,1])[k]
			te <- apply(tt,1,function(x) any(x==ll, na.rm=T))
			tef <- as.factor(te*1)
			levels(tef) <- c(levels(tef),"No sibilings or adopted")
			tef[rowSums(is.na(tt))==ncol(tt)] <- NA
			tef[bdE4$f.1767.0.0 == "Yes" | bdE4$f.1873.0.0b == 0 |  bdE4$f.1883.0.0b == 0] <- "No sibilings or adopted"
			bdE4[paste("f.20111.0.l",k,sep="")] <- tef
			levels(bdE4[[paste("f.20111.0.l",k,sep="")]])[2] <- ll
	}
bdE4 <- bdE4[,!(unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) =="20111" & unlist(lapply(strsplit(names(bdE4),"\\."),"[",4)) %in% as.character(seq(0,13,1)))]
# Exclude these variables because "Don't know", "Prefer not to answer", etc.. #
bdE4$f.20111.0.l1 <- NULL
bdE4$f.20111.0.l2 <- NULL
bdE4$f.20111.0.l3 <- NULL
bdE4$f.20111.0.l4 <- NULL
bdE4$f.20111.0.l5 <- NULL
bdE4$f.20111.0.l6 <- NULL

# EXCLUDE ILLNESS OF ADOPTED MOTHER, FATHER AND SIBILING - TOO FEW RECORDS #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("20112","20113","20114")]

# Country of Birth (non−UK origin) #
bdE4$f.20115.0.0 <- NULL

# Current employment status - corrected #
bdE4$f.20119.0.0 <- NULL

# BMI - duplicated, remove, other variable is better #
bdE4$f.21001.0.0 <- NULL

# Cascot confidence score - too many missing, no clear inclusion criteria #
bdE4$f.20121.0.0 <- NULL

# Weight (and manual weight) - duplicated, remove #
bdE4$f.21002.0.0 <- NULL
bdE4$f.3160.0.0 <- NULL

# Age when attended assessment centre - coded as special #
bdE4$f.21003.0.0 <- NULL

# Basophils count - Manual split in two categories or it wouldn't work #
bdE4$f.30160.0.0 <- as.factor(ifelse(bdE4$f.30160.0.0>5,"> 5","<= 5"))

# Nucleated red blood cell - Manual split in two categories or it wouldn't work #
bdE4$f.30170.0.0 <- as.factor(ifelse(bdE4$f.30170.0.0>5,"> 5","<= 5"))

# Nucleated red blood cell percentage - Manual split in two categories or it wouldn't work #
bdE4$f.30230.0.0 <- as.factor(ifelse(bdE4$f.30230.0.0>0.5,"> 0.5","<= 0.5"))

# Death date - set before as external variable #
bdE4$f.40000.0.0 <- NULL

# EXCLUDE ALL THESE REGISTRY VARIABLES #
bdE4 <- bdE4[,!unlist(lapply(strsplit(names(bdE4),"\\."),"[",2)) %in% 
c("40002","40005","40006","40008","40013","41078","41080","41083","41084","41085","41086","41087","41088","41089","41090","41091","41092","41093","41095","41096","41097","41098","41099","41100","41101","41102","41103","41104","41106","41107","41108","41109","41110","41111","41112","41132","41142","41146","41148","41142","41200","41201","41202","41204","41207","41208","41210","41219","41220","41221","41222","41223","41224","41225","41226","41227","41228","41229","41230","41231","41235")]

# EXCLUDE THESE VARIABLES BECAUSE FROM PILOT STUDY #
bdE4$f.40000.1.0 <- NULL
bdE4$f.40001.1.0 <- NULL

# N. of episodes of NHS episodes #
bdE4$f.41082.0.0 <- NULL

# N. of reported occurence of cancer #
bdE4$f.40009.0.0 <- NULL

# Age at death #
bdE4$f.40007.0.0 <- NULL

# Primary cause of death #
bdE4$f.40001.0.0 <- NULL

# REMOVER ALL f.xxxx.x.xb VARIABLES #
bdE4$f.3982.0.0b <- NULL
bdE4$f.3972.0.0b  <- NULL
bdE4$f.3983.0.0b  <- NULL
bdE4$f.3973.0.0b  <- NULL
bdE4$f.134.0.0b  <- NULL
bdE4$f.135.0.0b  <- NULL
bdE4$f.136.0.0b  <- NULL
bdE4$f.137.0.0b  <- NULL
bdE4$f.709.0.0b  <- NULL
bdE4$f.2734.0.0b  <- NULL
bdE4$f.3710.0.0b  <- NULL
bdE4$f.2804.0.0b  <- NULL
bdE4$f.1883.0.0b  <- NULL
bdE4$f.1873.0.0b  <- NULL
bdE4$f.1498.0.0b <- NULL
bdE4$f.1438.0.0b  <- NULL
bdE4$f.1458.0.0b  <- NULL

# Red Blood Cell Count - scale it because otherwise it creates computational problems #
bdE4$f.30010.0.0 <- scale(bdE4$f.30010.0.0)


# TRANSFORM ALL THE NUMERIC OR INTEGER CLASSES IN FACTORS #
for (i in 1:ncol(bdE4))
{
	if (colnames(bdE4)[i]%in% c("age","date_center","ddate","age","sex","center","f.eid","out","cofd","charlsonSR"))
	{
		bdE4[,i] <- bdE4[,i]
	}
	else if(class(bdE4[,i])=="ordered")
	{
		class(bdE4[,i]) <- "factor"
	}
	else if(class(bdE4[,i])=="numeric" | class(bdE4[,i])=="integer" | class(bdE4[,i])=="matrix")
	{
		bdE4[,i] <- cut2(as.numeric(bdE4[,i]))
		class(bdE4[,i]) <- "factor"
	}
	print(i)
}


## TRANSFORM EVERYTHING as FACTOR ##
class(bdE4$f.eid) <-  "numeric"
mode(bdE4$f.eid) <-  "numeric"
class(bdE4$age) <- "numeric"
mode(bdE4$age) <- "numeric"
class(bdE4$center) <- "factor"
class(bdE4$sex) <- "factor"
bdE4$charlsonSR <- as.factor(bdE4$charlsonSR)

## CHECK CLASS VARIABLES ##
#table(unlist(lapply(bdE4,class)))

# PLOT
set.seed(123)
ran <- sample(1:nrow(bdE4),10000)

pdf("Descriptives_statistics.pdf")
for (i in 1:ncol(bdE4))
{
	nameph <- ph[strsplit(names(bdE4)[i],"\\.")[[1]][2]==as.character(ph[,3]),2]
	if (sum(is.na(bdE4[,i]))<(nrow(bdE4)-30))
	{
		t <- bdE4[,i]
		if (class(t)=="Date" | class(t)=="numeric" | class(t)=="integer")
		{
			mat <- matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,0), 4, 4, byrow = TRUE)
			layout(mat)
			par(mar=c(2,2,2,2),oma = c(0, 0, 3, 0))
			hist(t, breaks=1000, ylab="",xlab="", main="")
			frame()
		legend("bottomleft",inset=0.3,c(paste("Min:",round(min(t,na.rm=T)),2),paste("Max:",round(max(t,na.rm=T)),2),paste("Mean:",round(mean(t, na.rm=T),2)),paste("Median:",round(median(t, na.rm=T),2)),paste("Missing:",sum(is.na(t))),paste("Non-Missing:",sum(!is.na(t))),paste("% Missing:",round(sum(is.na(t))/length(t)*100))))
			tryCatch(boxplot(t, ylab="",xlab=""),error = function(e) e)
			tryCatch(scatter.smooth(x=t[ran], y=bdE4$age[ran],xlab="", ylab="",pch=".", main="Age"),error = function(e) e)
			tryCatch(scatter.smooth(x=t[ran], y=bdE4$sex[ran],xlab="", ylab="",pch=".", main="Sex"),error = function(e) e)
			tryCatch(scatter.smooth(x=t[ran], y=bdE4$f.189.0.0[ran],xlab="", ylab="",pch=".", main="Deprivation"),error = function(e) e)
			mtext(paste(nameph,"-",names(bdE4)[i],"-",class(bdE4[,i])[1]), outer = TRUE)
		}
		if (class(t)=="ordered" | class(t)=="factor")
		{
			mat <- matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,0), 4, 4, byrow = TRUE)
			layout(mat)
			par(mar=c(2,4,2,2),oma = c(0, 0, 3, 0))
			xx <- barplot(table(t),ylab="",xlab="", main="", las=2)
			text(x = xx, y = 2000, label = round(table(t)/sum(table(t))*100,2), pos = 3, cex = 0.8, col = "red")
			text(x = xx, y = 15000, label = as.numeric(table(t)), pos = 3, cex = 0.8, col = "blue")
			frame()
	legend("bottomleft",inset=0.01,c(paste("Missing:",sum(is.na(t))),paste("Non-Missing:",sum(!is.na(t))),paste("% Missing:",round(sum(is.na(t))/length(t)*100))))
			frame()
			legend("topleft",col=c("blue","green","red"),c("1st tertile","2nd tertile","3rd tertile"), cex=0.8, pch=c(8,8,8))
			mtext(paste(nameph,"-",names(bdE4)[i],"-",class(bdE4[,i])[1]), outer = TRUE)
		}
	}
print(i)
}
dev.off()

bdE5 <- bdE4
# Save dataset #


save(bdE5,file="/proj/b2011036/uk.biobank/pre_imputation/out5.Rdata")