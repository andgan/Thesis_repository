/*----------------------------------------------
Filename: import_TWINGENE_25JUN2011.sas
Study: genscore
Author: Andrea Ganna
Date: 25JUN2011
Updated: 30JUN2011 - Duplicated MZ are now included
         11JUL2011 - Family history is included
         20SEP2011 - We keep 1 MZ per pair. I keep subject with chd=1 if pairs are discordant and 
                     randomly one of the two if both MZ have chd=0 or chd=1
		 11JUL2012 - Updated registries	and add outcomes
Purpose: Import data from twingene (ORACLE) for the genscore project, FHS variables and the GWAS key are included
Note:
-----------------------------------------------
Data used: databases from .twin
Data created: twge.twingene	 and twingene.txt
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/




**** ASSIGN LIBRARIES;

%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;
%let path_core=\\.psf\Home\Documents\Work\Phd_KI\Twingene\Data\Core_data;

*%let path=E:\Phd_KI\Genscore;
*%let path_core=E:\Phd_KI\Twingene\Data\Core_data;


* Settings libnames;
libname twge "&path.\Data\Twingene\Final";
libname core_tw "&path_core";

**********************************
**********************************
***           STEP 1           ***  
**********************************
********************************** 
*** CARDIOVASCULAR DisEASES    ***
***            AND             *** 
*** GENERAL INFO FROM TWINGENE ***
**********************************
**********************************


*** Enrollment date ***;
data date;
set core_tw.blood_sample;
check=mdy(substr(BLOOD_RECEIVED_DATE ,5,2),substr(BLOOD_RECEIVED_DATE ,7,2), substr(BLOOD_RECEIVED_DATE ,1,4));
keep TWINNR check;
format check date9.;
run;

*** Health questionnaire ***;

* Only subjects with checkdate;
data health;
set core_tw.health_questionnaire;
keep TWINNR BIRTHDATE DIABETES DIABETES_DOCTOR DIABETES_MEDICINE HEIGHT_CM 
     SYST_ONE_MMHG SYST_TWO_MMHG WEIGHT_KG HART_MEDICINE_TYPE;
run;


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
proc sort data=date; by twinnr; run;
data cvd;
merge cvd date (keep=twinnr check);
by twinnr;
if indat =. or check=. then delete;
run;


*** Cause of deaths ***;

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
proc sort data=date; by twinnr; run;
data death;
merge death date (keep=twinnr check);
by twinnr;
if deathd =. or check=. then delete;
if deathd => check ;
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


*** Diabetes from registry  ***;

** Drug registry;

data diab_drug;
set core_tw.Link_drug_prescription_info;
indat= datepart(FDATUM);
diab_reg=1;
where find (ATC, "A10")=1;
keep twinnr diab_reg indat;
format indat date9.;
run;

** Patient registry;

data diab;
set core_tw.v_link_patient_2012_andgan;
indat=mdy(substr(INDATE,5,2),substr(INDATE,7,2), substr(INDATE,1,4));
if indat=. then delete;
keep TWINNR indat HDIA ICD_VERSION;
format indat date9.;
run;

*** Coding ***;
data diab;
set diab;

if ICD_VERSION=10      and find(HDIA, "E11")=1 then diab_reg=1;
else if ICD_VERSION=10 and find(HDIA, "E14")=1 then diab_reg=1 ;
else if ICD_VERSION=9 and find(HDIA, "250")=1 then diab_reg=1;
else if ICD_VERSION=8 and find(HDIA, "250")=1 then diab_reg=1;
else if ICD_VERSION=7 and find(HDIA, "260")=1 then diab_reg=1;


if diab_reg=1;
drop  ICD_VERSION HDIA;
run;


data diab_m;
set diab_drug diab;
run;



*** Subjects with diabetes at baseline ***;

proc sort data=diab_m; by TWINNR indat; run;
proc sort data=date; by twinnr; run;
data diab_m;
merge diab_m date (keep=twinnr check);
by twinnr;
diabd = indat;
if check =. or diabd=. then delete;	
format diabd date9.;
run;


** Keep only first diabetes event;

proc sort data=diab_m; by twinnr diabd; run;
data diab_m;
set diab_m;
by twinnr;
if first.twinnr;
keep twinnr diab_reg diabd;
run;


*** Sex and PAIRID***;

data sex;
set core_tw.V_admin;
keep TWINNR SEX PAIRID;
run;

*** Laboratory data ***;

data labdata;
set core_tw.V_labdata;
keep TWINNR ANALYSisDESC VALUE;
run;
proc sort data=labdata; by twinnr; run; 
proc transpose data=labdata out=labdata ; by twinnr; var value; id ANALYSisDESC; run;
data labdata; set labdata; keep TWINNR S_Kolesterol S_HDL_Kolesterol B_HbA1c; run;

****** DATA FROM SALT *******;

data smoke;
set core_tw.X_ibs_az;
keep twinnr XCSMOKER;
run;


***** GWAS KEY ****;

* Not considering the duplicated MZ;
data key;
set core_tw.gwas_key;
*if IMPUTE_MZ=0 ** CHANGED 30JUN2011 **;
drop IMPUTE_MZ PLINK_KEY;
run;

**********************************
**********************************
***           STEP 4           ***  
**********************************
********************************** 
***       GREAT MERGING        ***
**********************************
**********************************;


****** MERGE *****;
proc sort data=health; by twinnr; run;
proc sort data=cad; by twinnr; run;  
proc sort data=chd; by twinnr; run; 
proc sort data=is; by twinnr; run; 
proc sort data=hs; by twinnr; run; 
proc sort data=hf; by twinnr; run; 
proc sort data=death; by twinnr; run; 
proc sort data=diab_m; by twinnr; run; 
proc sort data=sex; by twinnr; run; 
proc sort data=labdata; by twinnr; run; 
proc sort data=smoke; by twinnr; run;
proc sort data=key; by twinnr; run;
proc sort data=date; by twinnr; run;




data twingene;
merge health cad chd is hs hf death diab_m sex labdata smoke date key (in=x);
by twinnr;
if x=1;
run;

*** Now for the MZ with two pairs I selected the subject with the CHD event ;

* Keep onlt MZ with 2 pairs;
proc freq data=twingene noprint;
table  GWAS_FID / out=temp;
where bestzyg=1;
run;
proc sort data=twingene; by GWAS_FID;
data mz_1;
merge twingene temp;
by GWAS_FID;
if COUNT=2;
run;

* Exclude=1 if pairs are discordant and chd=0 or randomly one of the two if both chd=0 or chd=1;
proc sort data=mz_1; by GWAS_FID chd;
data MZ;
set mz_1;
by gwas_fid;
if first.GWAS_FID then exclude=1;
run;

* Remerge with twingene;
proc sort data=twingene; by GWAS_ID;
proc sort data=MZ; by GWAS_ID;
data twingene;
merge twingene MZ;
by GWAS_ID;
if exclude ne 1;
drop exclude PERCENT COUNT;
run;

*************************************
***** CREATING TWINGENE DATASET *****
*************************************;

data twge.twingene;
set twingene;

* From alphanumeric to numeric;
diabetess=input(DIABETES, best12.);
diabetess_doctor=input(DIABETES_DOCTOR, best12.);
diabetess_medicine=input(DIABETES_MEDICINE, best12.);
height=input(HEIGHT_CM, best12.);
SYST_ONE_MMHG_=input(SYST_ONE_MMHG, best12.);
SYST_TWO_MMHG_=input(SYST_TWO_MMHG, best12.);
weight=input(WEIGHT_KG, best12.);
tc=input(S_Kolesterol, best12.);
hdl=input(S_HDL_Kolesterol, best12.);
HbA1c=input(B_HbA1c, best12.);

* From alphanumeric to date;
birth_d=mdy(substr(BIRTHDATE,5,2),substr(BIRTHDATE,7,2), substr(BIRTHDATE,1,4));
check_d=check;
death_d=deathd;


* Missing;

if weight in (999,998,997) then weight=.;
if height in (999,998,997) then height =.;
if diabetess in (999,998,997) then diabetess =.;
if diab_reg =. then diab_reg=0;
if cad =. then cad=0;
if chd =. then chd=0;
if is =. then is=0;
if hs =. then hs=0;
if hf =. then hf=0;
if death =. then death=0;


* Medications;

* Beta blockers;
if index (HART_MEDICINE_TYPE, "Tenormin")>0 then betab=1;
else if index (HART_MEDICINE_TYPE, "Atenolol")>0 then betab=1;
else if index (HART_MEDICINE_TYPE, "Seloken")>0 then betab=1;
else betab=0;
*ACE inibitors/Angiotensin II-inhibitors;
if index (HART_MEDICINE_TYPE, "Pramace")>0 then ACE_angio=1;
else if index (HART_MEDICINE_TYPE, "Renitec")>0 then ACE_angio=1;
else if index (HART_MEDICINE_TYPE, "Triatec")>0 then ACE_angio=1;
else if index (HART_MEDICINE_TYPE, "Cozaar")>0 then ACE_angio=1;
else if index (HART_MEDICINE_TYPE, "Enalapril")>0 then ACE_angio=1;
else if index (HART_MEDICINE_TYPE, "Zestril")>0 then ACE_angio=1;
else ACE_angio=0;
*Calcium blockers;
if index (HART_MEDICINE_TYPE, "Cardizem")>0 then calcb=1;
else if index (HART_MEDICINE_TYPE, "Plendil")>0 then calcb=1;
else calcb=0;
*Diuretics;
if index (HART_MEDICINE_TYPE, "Lasix")>0 then diur=1;
else if index (HART_MEDICINE_TYPE, "Salures")>0 then diur=1;
else diur=0;
*Statins;
if index (HART_MEDICINE_TYPE, "Lipitor")>0 then statins=1;
else if index (HART_MEDICINE_TYPE, "Simvastatin")>0 then statins=1;
else if index (HART_MEDICINE_TYPE, "Zocord")>0 then statins=1;
else statins=0;
*Fibrates;
if index (HART_MEDICINE_TYPE, "Bezalip")>0 then fibr=1;
else fibr=0;

* Define Antihypertensive treatment;
if betab=1 or ACE_angio=1 or calcb=1 or diur=1 then antihyp=1;
else antihyp=0;

* Define lowering lipid drougs;
if fibr=1 or statins=1 then drug_lipid=1;
else drug_lipid=0;

* Only real checkdate;
if check_d>18627 or check_d=. then delete; * 31DEC2010;


* Setting variables;
if SYST_ONE_MMHG_ ne . and SYST_TWO_MMHG_ ne . then
 sbp=(SYST_ONE_MMHG_+SYST_TWO_MMHG_)/2;
else if SYST_ONE_MMHG_ ne . and SYST_TWO_MMHG_ eq . then 
 sbp=(SYST_ONE_MMHG_);


smoke=XCSMOKER;

if weight ne . and height ne . then
 bmi=weight/(height/100)**2;
else bmi=.;

age=(check_d-birth_d)/365.25;


*** DEFINE OUTCOMES;

* Define diabetes (at baseline);

if (diabetess=1 and diabetess_doctor=1) or diabetess_medicine=1 then diab=1; else diab=0;


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


study="twingene";

keep GWAS_ID GWAS_FID bestzyg check_d birth_d age sex hdl tc bmi sbp smoke antihyp diab cad chd is hs hf inccad incis incchd inchs inchf 
cadd chdd isd hsd hfd age_entry_cad age_exit_cad age_entry_chd age_exit_chd age_entry_is age_exit_is age_entry_hs age_exit_hs age_entry_hf age_exit_hf
TWINNR study drug_lipid;

format birth_d date9. check_d date9. chdd date9. isd date9.; 
run;


** Calculate family history as:
1 both twins had a CHD or 1 twin has a chd and the other not (then the twin without chd has famhist=1 and the other famhist=0)
0 if both twins does not have a chd
. if only 1  twin is in twingene.;

proc sort data=twge.twingene; by GWAS_FID; run;
proc means data=twge.twingene noprint max;
by GWAS_FID;
var chd;
output out=temp sum=sum;
run;

data twge.twingene;
merge twge.twingene temp;
by GWAS_FID;
if _FREQ_=2 and sum=2 then 	famhist=1;
else if _FREQ_=2 and sum=1 and chd=0 then famhist=1;
else if _FREQ_=2 and sum=1 and chd=1 then famhist=0;
else if _FREQ_=2 and sum=0 then famhist=0;
else famhist=.;

drop _TYPE_ _FREQ_ sum;
run;


PROC EXPORT DATA= twge.twingene
            OUTFILE= "&path.\Data\Twingene\Final\twingene.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=YES;
RUN;


**Descriptive;


proc freq data=twge.twingene;
table smoke antihyp diab ;
run;

proc means data=twge.twingene;
var bmi sbp tc hdl;
run;


** Association;

proc phreg data=twge.twingene;
model (age_entry_chd,age_exit_chd)*incchd(0) = sex tc hdl antihyp sbp smoke bmi diab famhist;
where chdd > check_d;
run;

