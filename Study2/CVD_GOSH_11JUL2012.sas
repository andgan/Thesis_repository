/*----------------------------------------------
Filename: CVD_GOSH.sas
Study: genscore
Author: Andrea Ganna
Date: 29JUN2011
Updated: 11JUL2011 - Included family history informations
		 11JUL2012 - Updated the registries	and add outcomes
Purpose: We attach to GOSH data CVD outcomes form the registries;
Note:
-----------------------------------------------
Data used: gosh.GOSH core_tw data from the registries
Data created: gosh.GOSH_CVD
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/




**** ASSIGN LIBRARIES;

%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;
%let path_core=\\.psf\Home\Documents\Work\Phd_KI\Twingene\Data\Core_data;

*%let path=E:\Phd_KI\Genscore;
*%let path_core=E:\Phd_KI\Twingene\Data\Core_data;


* Settings libnames;
libname gosh "&path.\Data\GOSH\Final";
libname core_tw "&path_core";



*** Incident CVD events ***;

data cvd;
set core_tw.v_link_patient_2012_andgan;
indat=mdy(substr(INDATE,5,2),substr(INDATE,7,2), substr(INDATE,1,4));
if indat=. then delete;
keep TWINNR indat HDIA OP1-OP12 ICD_VERSION;
format indat date9.;
run;


*** Coding ***;
data cvd;
set cvd;
array op_a {12} OP1-OP12;

** CAD (MI + surgical codes);

if ICD_VERSION=10 and find(HDIA, "I21")=1 then hcad=1 ;
else if ICD_VERSION=10 and find(HDIA, "I22")=1 then hcad=1;
else if ICD_VERSION=9 and find(HDIA, "410")=1 then hcad=1;
else if ICD_VERSION=8 and find(HDIA, "410")=1 then hcad=1;


do i=1 to 12;

if op_a {i} in ("FNG02","FNG05") then hcad=1;
else if  find(op_a {i}, "FNC")=1 then hcad=1;
else if  find(op_a {i}, "FND")=1 then hcad=1;
else if  find(op_a {i}, "FNE")=1 then hcad=1;
else if  find(op_a {i}, "3080")=1 then hcad=1;
else if  find(op_a {i}, "3127")=1 then hcad=1;
else if  find(op_a {i}, "3158")=1 then hcad=1;
end;
if hcad =. then hcad=0;


** CHD;

if ICD_VERSION=10      and find(HDIA, "I200")=1 then hchd=1;
else if ICD_VERSION=10 and find(HDIA, "I21")=1 then hchd=1 ;
else if ICD_VERSION=10 and find(HDIA, "I22")=1 then hchd=1;
else if ICD_VERSION=9 and find(HDIA, "410")=1 then hchd=1;
else if ICD_VERSION=9 and find(HDIA, "411B")=1 then hchd=1;
else if ICD_VERSION=8 and find(HDIA, "410")=1 then hchd=1;
else if ICD_VERSION=8 and find(HDIA, "411")=1 then hchd=1;


do i=1 to 12;
if op_a {i} in ("FNG02","FNG05") then hchd=1;
else if  find(op_a {i}, "FNC")=1 then hchd=1;
else if  find(op_a {i}, "FND")=1 then hchd=1;
else if  find(op_a {i}, "FNE")=1 then hchd=1;
else if  find(op_a {i}, "3080")=1 then hchd=1;
else if  find(op_a {i}, "3127")=1 then hchd=1;
else if  find(op_a {i}, "3158")=1 then hchd=1;
end;
if hchd =. then hchd=0;



** Ischemic Stroke;

if ICD_VERSION=10     and find(HDIA, "I63")=1 then his=1;
else if ICD_VERSION=9 and find(HDIA, "433")=1 then his=1 ;
else if ICD_VERSION=9 and find(HDIA, "434")=1 then his=1;
else if ICD_VERSION=8 and find(HDIA, "432")=1 then his=1;
else if ICD_VERSION=8 and find(HDIA, "433")=1 then his=1;
else if ICD_VERSION=8 and find(HDIA, "434")=1 then his=1;
else his=0;



** Hemorrhagic Stroke;

if ICD_VERSION=10     and find(HDIA, "I60")=1 then hhs=1;
else if ICD_VERSION=10 and find(HDIA, "I61")=1 then hhs=1;
else if ICD_VERSION=10 and find(HDIA, "I62")=1 then hhs=1;
else if ICD_VERSION=9 and find(HDIA, "430")=1 then hhs=1 ;
else if ICD_VERSION=9 and find(HDIA, "431")=1 then hhs=1;
else if ICD_VERSION=9 and find(HDIA, "432")=1 then hhs=1;
else if ICD_VERSION=8 and find(HDIA, "430")=1 then hhs=1;
else if ICD_VERSION=8 and find(HDIA, "431")=1 then hhs=1;
else hhs=0;



** Heart Failure;

if ICD_VERSION=10     and find(HDIA, "I50")=1 then hhf=1;
else if ICD_VERSION=9 and find(HDIA, "428")=1 then hhf=1 ;
else if ICD_VERSION=8 and find(HDIA, "427.00")=1 then hhf=1;
else if ICD_VERSION=8 and find(HDIA, "427.10")=1 then hhf=1;
else hhf=0;


if hchd=1 or his=1 or hhs=1 or hhf=1 or hcad=1;


drop  ICD_VERSION i;
run;

* Select subjects in TWINGENE;

proc sort data=cvd; by TWINNR; run;
proc sort data=gosh.GOSH; by twinnr; run;
data cvd;
merge cvd gosh.GOSH (keep=twinnr in=x);
by twinnr;
if indat =. then delete;
if x=1;
run;


*** Cuase of deaths ***;

data cdeath;
set core_tw.link_death_cause_2012;
if MORSAK_NR=1 then output;
run;

*** Deaths ***;

data adeath;
set core_tw.V_person_complete;
deathd=mdy(substr(DEATHDATE ,5,2),substr(DEATHDATE ,7,2), substr(DEATHDATE ,1,4));
keep TWINNR deathd;
if deathd ne . and deathd <= 18627; * 31DEC2010;
format deathd date9.;
run;


proc sort data=adeath; by twinnr; run;
proc sort data=cdeath; by twinnr; run;

data death;
merge adeath cdeath;
by twinnr;

** CAD;

if find(MORSAK, "I21")=1 then dcad=1 ;
else if find(MORSAK, "I22")=1 then dcad=1;
else if find(MORSAK, "410")=1 then dcad=1;
else dcad=0;



** CHD;

if find (MORSAK ,"I200") then dchd=1;
else if find (MORSAK, "411B") then dchd=1;
else if find(MORSAK, "I21")=1 then dchd=1 ;
else if find(MORSAK, "I22")=1 then dchd=1;
else if find(MORSAK, "410")=1 then dchd=1;
else dchd=0;


** Ischemic Stroke;

if find(MORSAK, "I63")=1 then dis=1;
else if find(MORSAK, "433")=1 then dis=1;
else if find(MORSAK, "434")=1 then dis=1;
else dis=0;


** Hemorrhagic Stroke;


if find(MORSAK, "I60")=1 then dhs=1;
else if find(MORSAK, "I61")=1 then dhs=1;
else if find(MORSAK, "I62")=1 then dhs=1;
else if find(MORSAK, "430")=1 then dhs=1 ;
else if find(MORSAK, "431")=1 then dhs=1;
else if find(MORSAK, "432")=1 then dhs=1;
else dhs=0;


** Heart Failure;

if find(MORSAK, "I50")=1 then dhf=1;
else if find(MORSAK, "428")=1 then dhf=1 ;
else if find(MORSAK, "427.00")=1 then dhf=1;
else if find(MORSAK, "427.10")=1 then dhf=1;
else dhf=0;



* deleting subject with cause of death but not death date;
if deathd =. then delete;
if dchd ne 1 and dis ne 1 and dhs ne 1 and dhf ne 1 and dcad ne 1 then dother=1;
else dother=0;

drop MORSAK_NR;
run;


* Check no patients checked after death;

proc sort data=death; by TWINNR; run;
proc sort data=gosh.GOSH; by twinnr; run;
data death;
merge death gosh.GOSH (keep=twinnr check_d in=x);
by twinnr;
if x=1;
if deathd => check_d and deathd ne .;
run;


 *** Merge hospital registry and death registries;
proc sort data=death; by twinnr; run;
proc sort data=cvd; by twinnr; run;

data cvddea;
merge cvd death;
by twinnr;

* Define main outcome;
if hcad =1 or dcad =1 then cad = 1; else cad=0;
if hchd =1 or dchd =1 then chd = 1; else chd=0;
if his = 1 or dis = 1 then is = 1; else is=0;
if hhs = 1 or dhs = 1 then hs = 1; else hs=0;
if hhf = 1 or dhf = 1 then hf = 1; else hf=0;
if deathd ne . then death=1; else death=0;


* Define events dates;
if hcad=1  then cadd=indat;
else if dcad=1  then cadd=deathd;
if hchd=1  then chdd=indat;
else if dchd=1  then chdd=deathd;
if his = 1  then isd=indat;
else if dis=1  then isd=deathd;
if hhs = 1  then hsd=indat;
else if dhs=1  then hsd=deathd;
if hhf = 1  then hfd=indat;
else if dhf=1  then hfd=deathd;


format chdd date9. isd date9. hsd date9. hfd date9. cadd date9.; 

run;



*** GET FIRST CAD date;

proc sort data=cvddea; by twinnr cadd; run;

data cad;
set cvddea;
by twinnr;
if first.twinnr;
keep twinnr cad cadd;
where cad=1;
run;


*** GET FIRST CHD date;

proc sort data=cvddea; by twinnr chdd; run;


data chd;
set cvddea;
by twinnr;
if first.twinnr;
keep twinnr chd chdd;
where chd=1;
run;


*** GET FIRST IS date;

proc sort data=cvddea; by twinnr isd; run;


data is;
set cvddea;
by twinnr;
if first.twinnr;
keep twinnr is isd;
where is=1;
run;



*** GET FIRST HS date;

proc sort data=cvddea; by twinnr hsd; run;


data hs;
set cvddea;
by twinnr;
if first.twinnr;
keep twinnr hs hsd;
where hs=1;
run;



*** GET FIRST HF date;

proc sort data=cvddea; by twinnr hfd; run;


data hf;
set cvddea;
by twinnr;
if first.twinnr;
keep twinnr hf hfd;
where hf=1;
run;



* Get death date;

proc sort data=cvddea; by twinnr deathd; run;

data death;
set cvddea;
by twinnr;
if  first.twinnr;
keep twinnr death deathd;
where death=1;
run;


proc sort data=cad; by twinnr; run;  
proc sort data=chd; by twinnr; run; 
proc sort data=is; by twinnr; run; 
proc sort data=hs; by twinnr; run; 
proc sort data=hf; by twinnr; run; 
proc sort data=death; by twinnr; run; 
proc sort data=gosh.GOSH; by twinnr; run;




data gosh_CVD;
merge cad chd is hs hf death gosh.GOSH (in=x);
by twinnr;
if x=1;
run;


*************************************
***** CREATING GOSH_CVD DATASET *****
*************************************;

data gosh.GOSH_CVD;
set gosh_CVD;

if cad =. then cad=0;
if chd =. then chd=0;
if is =. then is=0;
if hs =. then hs=0;
if hf =. then hf=0;
if death =. then death=0;

age=(check_d-birth_d)/365.25;


** Event date;
if cad ne 1 and death=1 then cadd=deathd;
else if cad ne 1 and death ne 1 then cadd=18627;

if chd ne 1 and death=1 then chdd=deathd;
else if chd ne 1 and death ne 1 then chdd=18627;

if is ne 1 and death=1 then isd=deathd;
else if is ne 1 and death ne 1 then isd=18627;

if hs ne 1 and death=1 then hsd=deathd;
else if hs ne 1 and death ne 1 then hsd=18627;

if hf ne 1 and death=1 then hfd=deathd;
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

* Creat these variable to be compatible with the other studies;
GWAS_ID = input(TWINNR, $20.);
GWAS_FID = input(PAIRID, $20.);

*study="GOSH";

keep GWAS_ID GWAS_FID bestzyg check_d birth_d age sex hdl tc bmi sbp smoke antihyp diab cad chd is hs hf inccad incis incchd inchs inchf 
cadd chdd isd hsd hfd age_entry_cad age_exit_cad age_entry_chd age_exit_chd age_entry_is age_exit_is age_entry_hs age_exit_hs age_entry_hf age_exit_hf
study drug_lipid;

format  cadd date9. chdd date9. isd date9. hsd date9. hfd date9.; 
run;

** Calculate family history as:
1 both twins had a CHD or 1 twin has a chd and the other not (then the twin without chd has famhist=1 and the other famhist=0)
0 if both twins does not have a chd
. if only 1  twin is in twingene.;

proc sort data=gosh.GOSH_CVD; by GWAS_FID; run;
proc means data=gosh.GOSH_CVD noprint max;
by GWAS_FID;
var chd;
output out=temp sum=sum;
run;

data gosh.GOSH_CVD;
merge gosh.GOSH_CVD temp;
by GWAS_FID;
if _FREQ_=2 and sum=2 then 	famhist=1;
else if _FREQ_=2 and sum=1 and chd=0 then famhist=1;
else if _FREQ_=2 and sum=1 and chd=1 then famhist=0;
else if _FREQ_=2 and sum=0 then famhist=0;
else famhist=.;

drop _TYPE_ _FREQ_ sum;
run;




PROC EXPORT DATA= gosh.GOSH_CVD
            OUTFILE= "&path.\Data\Twingene\Final\gosh_CVD.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=YES;
RUN;


** check is association makes sense;

proc phreg data=gosh.GOSH_CVD;
model (age_entry_chd,age_exit_chd)*incchd(0) = sex tc hdl antihyp sbp smoke bmi diab famhist;
where chdd > check_d;
run;
