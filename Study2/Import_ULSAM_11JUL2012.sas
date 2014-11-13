/*----------------------------------------------
Filename: import_ULSAM.sas
Study: genscore
Author: Andrea Ganna
Date: 25JUN2011
Updated: 01JUL2011 - Now pat is called GWAS_ID and is done as: "UL"||pat
         04JUL2011 - Now added the .fam file, so only genotyped subjects are included
         20SEP2011 - Added lipid lowering drugs
         20OCT2011 - Modify Ulsam_70
		 11JUL2012 - Updated registries and add outcomes
Purpose: Import data from ULSAM for the genscore project. Create ULSAM_50 and ULSAM_70
Note:
-----------------------------------------------
Data used: hosp110623.dta death110623.dta d110621.dta ingel70+register2002-modified.dta
Data created: ULSAM_50.txt ULSAM_70.txt
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/



%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;

libname ulsam "&path.\Data\ULSAM\Final";

PROC IMPORT OUT= WORK.hosp 
            DATAFILE= "&path.\Data\ULSAM\Original\hosp111013.dta" 
            DBMS=DTA REPLACE;
RUN;

PROC IMPORT OUT= WORK.death 
            DATAFILE= "&path.\Data\ULSAM\Original\death111013.dta" 
            DBMS=DTA REPLACE;
RUN;

PROC IMPORT OUT= WORK.med_50 
            DATAFILE= "&path.\Data\ULSAM\Original\allvar50+allHDR-validerad-loggad-extravariabler-2001.dta" 
            DBMS=DTA REPLACE;
RUN;

* Keep only medications;
data med_50; set med_50; keep pat x100 x303; run;

*** Merge hospital registry with death registry and keep the first event date and death date;

proc sort data=hosp; by pat; run;
proc sort data=death; by pat; run;

data hosdea;
merge hosp death;
by pat;
evdate = mdy(substr(INDATE ,5,2),substr(INDATE ,7,2), substr(INDATE ,1,4));
ddate =  mdy(substr(P013 ,5,2),substr(P013 ,7,2), substr(P013 ,1,4));

* Death after 31dec2010 are not considered;
if ddate > 18627 then ddate = .;

* Define main outcome;


if CHD_HOSP in (1,3) then cad_hosp = 1; else cad_hosp=0;
if (ICDVER=10 and find(HDIA, "I200")=1) or (ICDVER=9 and find(HDIA, "411B")=1) or (ICDVER=8 and find(HDIA, "411")=1) then cad_hosp=0;


if CHD_DEATH =1 then cad_death = 1; else cad_death=0;
if (ICDVER=10 and find(P002, "I200")=1) or (ICDVER=9 and find(P002, "411B")=1) or (ICDVER=8 and find(P002, "411")=1) then cad_hosp=0;

if cad_hosp in (1,3) or cad_death =1 then cad = 1; else cad=0;
if CHD_HOSP in (1,3) or CHD_DEATH =1 then chd = 1; else chd=0;
if ISTROKE_HOSP = 1 or ISTROKE_DEATH = 1 then is = 1; else is=0;
if HSTROKE_HOSP = 1 or HSTROKE_DEATH = 1 then hs = 1; else hs=0;
if HFAIL_HOSP = 1 or HFAIL_DEATH = 1 then hf = 1; else hf=0;

if ddate ne . then death=1; else death=0;

* Define events dates;
if  cad_hosp=1 then cadd=evdate;
else if cad_death=1  then cadd=ddate;

if  CHD_HOSP in (1,3) then chdd=evdate;
else if CHD_DEATH=1  then chdd=ddate;

if ISTROKE_HOSP = 1 then isd=evdate;
else if ISTROKE_DEATH=1  then isd=ddate;

if HSTROKE_HOSP = 1 then hsd=evdate;
else if HSTROKE_DEATH=1  then hsd=ddate;

if HFAIL_HOSP = 1 then hfd=evdate;
else if HFAIL_DEATH=1  then hfd=ddate;


format evdate date9. ddate date9. cadd date9. chdd date9. isd date9. hsd date9. hfd date9.; 
run;



*** GET FIRST CAD date;
proc sort data=hosdea; by pat cadd; run;

data cad;
set hosdea;
by pat;
if first.pat;
keep pat cad cadd;
where cad=1;
run;

*** GET FIRST CHD date;
proc sort data=hosdea; by pat chdd; run;

data chd;
set hosdea;
by pat;
if first.pat;
keep pat chd chdd;
where chd=1;
run;


*** GET FIRST IS date;
proc sort data=hosdea; by pat isd; run;

data is;
set hosdea;
by pat;
if first.pat;
keep pat is isd;
where is=1;
run;

*** GET FIRST HS date;
proc sort data=hosdea; by pat hsd; run;

data hs;
set hosdea;
by pat;
if first.pat;
keep pat hs hsd;
where hs=1;
run;

*** GET FIRST HF date;
proc sort data=hosdea; by pat hfd; run;

data hf;
set hosdea;
by pat;
if first.pat;
keep pat hf hfd;
where hf=1;
run;


* Get death date;
proc sort data=hosdea; by pat ddate; run;
data death;
set hosdea;
by pat;
if  first.pat;
keep pat death ddate;
where death=1;
run;



**IMPORT fam file so i keep only these subjects;
PROC IMPORT OUT= WORK.Fam 
            DATAFILE= "&path.\Data\ULSAM\ulsam.+.qced.imputed.ceu.p2.fam" 
            DBMS=DLM REPLACE;
     DELIMITER=' '; 
     GETNAMES=NO;
RUN;



** Import phenotype informations;

PROC IMPORT OUT= WORK.ulsam_other 
            DATAFILE= "&path.\Data\ULSAM\Original\d110621.dta" 
            DBMS=DTA REPLACE;
RUN;


PROC IMPORT OUT= WORK.ULSAM_sel 
            DATAFILE= "&path.\Data\ULSAM\Original\ingel70+register2002-modified.dta" 
            DBMS=DTA REPLACE;
RUN;



proc sort data=ulsam_other; by pat; run;
proc sort data=ulsam_sel; by pat; run;
proc sort data=fam; by VAR2; run;
proc sort data=med_50; by pat; run;

*** MERGING ULSAM_70;

data ulsam.ulsam_70;
merge cad chd is hs hf death ulsam_other ( keep=pat Z116-Z135 DAT70) ulsam_sel (rename=(diab=diab_)) fam (rename=(VAR2=pat) in=x);
by pat;

* Rename;
diab_=Z378;
diab_treat=z405;
glycemia=z319;
hdl=z302;
tc=z972;
age=age70;
sbp=z013;
antihyp=z101;
smoke= z085;
bmi=z290;


**Dates;
bdate_=input(bdate, $12.);
diabdat_=input(diabdat, $12.);

check_d = mdy(substr(DAT70 ,5,2),substr(DAT70 ,7,2), substr(DAT70 ,1,4)); 
birth_d =  mdy(substr(bdate_ ,5,2),substr(bdate_ ,7,2), substr(bdate_ ,1,4));
diab_d =  mdy(substr(diabdat_ ,5,2),substr(diabdat_ ,7,2), substr(diabdat_ ,1,4));


**Diabetes;
if diab_ = . and diab_treat=. then diab=.;
else if (diab_ = 1 and diab_d < check_d) or diab_treat = 1  then diab=1;
else diab=0;

** Lipid lowering drugs;
if Z105=1 then drug_lipid=1; 
else drug_lipid=0;

** Event date;
if cad ne 1 and death=1 then cadd=ddate;
else if cad ne 1 and death ne 1 then cadd=18627;

if chd ne 1 and death=1 then chdd=ddate;
else if chd ne 1 and death ne 1 then chdd=18627;

if is ne 1 and death=1 then isd=ddate;
else if is ne 1 and death ne 1 then isd=18627;

if hs ne 1 and death=1 then hsd=ddate;
else if hs ne 1 and death ne 1 then hsd=18627;

if hf ne 1 and death=1 then hfd=ddate;
else if hf ne 1 and death ne 1 then hfd=18627;


** Incident events;
if cadd > check_d and cad=1 then inccad=1; else inccad=0;
if chdd > check_d and chd=1 then incchd=1; else incchd=0;
if isd > check_d and is=1 then incis=1; else incis=0;
if hsd > check_d and hs=1 then inchs=1; else inchs=0;
if hfd > check_d and hf=1 then inchf=1; else inchf=0;


** Age entry and exit;
age_entry_cad = age;
age_exit_cad = 	(cadd - birth_d) / 365.25;

age_entry_chd = age;
age_exit_chd = 	(chdd - birth_d) / 365.25;

age_entry_is = age;
age_exit_is = (isd - birth_d) / 365.25;

age_entry_hs = age;
age_exit_hs = (hsd - birth_d) / 365.25;

age_entry_hf = age;
age_exit_hf = (hfd - birth_d) / 365.25;

* outcomes;
if cad =. then cad=0;
if chd =. then chd=0;
if is =. then is=0;
if hs =. then hs=0;
if hf =. then hf=0;

* Family history;

if Z117=. and Z118=. and  Z119=. and  Z120=. and  Z121=. and  Z122=. and  Z123=. and  Z124=. and 
Z125=. and  Z126=. and  Z127=. and  Z132=. and  Z133=. and  Z134=. and  Z135=. then  famhist=.;
else if Z117-Z127 = 1 or Z132-Z135 =1 then famhist=1;
else famhist=0;

* Sex;
sex=1;

* Keep only subjects with checkdate and only with genotype information;
if check_d ne . ;
if x=1;

* To make it comparable;
GWAS_ID="UL"||input(pat,$20.);
GWAS_FID="UL"||input(pat,$20.);

study="ulsam";

keep GWAS_ID GWAS_FID check_d birth_d age sex hdl tc bmi sbp smoke antihyp famhist diab cad chd is hs hf inccad incis incchd inchs inchf 
cadd chdd isd hsd hfd age_entry_cad age_exit_cad age_entry_chd age_exit_chd age_entry_is age_exit_is age_entry_hs age_exit_hs age_entry_hf age_exit_hf
study drug_lipid;

format check_d date9. birth_d date9.;
run;



PROC EXPORT DATA= ULSAM.ULSAM_70
            OUTFILE= "&path.\Data\ULSAM\Final\ULSAM_70.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=YES;
RUN;


* Ulsam 70;

proc freq data=ulsam.ulsam_70;
table smoke antihyp diab ;
run;

proc means data=ulsam.ulsam_70;
var bmi sbp tc hdl;
run;


** Association;

proc phreg data=ulsam.ulsam_50;
model (age_entry_chd,age_exit_chd)*incchd(0) = sex tc hdl antihyp sbp smoke bmi diab famhist;
where chdd > check_d;
run;
