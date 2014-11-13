/*----------------------------------------------
Filename: import_TWINGENE.sas
Study: CCNCC
Author: Andrea Ganna
Date: 15JUL2010
Updated: 30SEP2010
Purpose: Import data from twingene (ORACLE) and export data for R
Note:
-----------------------------------------------
Data used: ORACLE databases 
Data created: twingene.sas7bdat
              twingene_unr.sas7bdat
              wc.csv (for R use)
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/



*Setting paths;

%let path=\\.psf\Home\Documents\Work\Phd_KI\CCNCC;
*%let path=E:\Phd_KI\CCNCC;

* Settings libnames;

libname tw "&path.\Data";



****** DATA FROM TWINGENE ******;


*** Health questionnaire ***;

* Only subjects with checkdate;
data health;
set twin.health_questionnaire;
check=mdy(substr(CHECKDATE ,5,2),substr(CHECKDATE ,7,2), substr(CHECKDATE ,1,4));
if check = . then delete; *Deleting twins without checkdate;
keep TWINNR BIRTHDATE DIABETES DIAS_ONE_MMHG DIAS_TWO_MMHG HEIGHT_CM 
     SYST_ONE_MMHG SYST_TWO_MMHG WEIGHT_KG HART_MEDICINE_TYPE check;
format check date9.;
run;


*** Incident CVD events ***;

data cvd;
set str.link_patient_2010;
indat=mdy(substr(INDATE,5,2),substr(INDATE,7,2), substr(INDATE,1,4));
if indat=. then delete;
keep TWINNR indat HDIA OP1-OP3 ICD_VERSION;
format indat date9.;
run;


*** Coding ***;
data cvd;
set cvd;
array op_a {3} OP1-OP3;

if ICD_VERSION=10      and find(HDIA, "I200")=1 then CHD=1;
else if ICD_VERSION=10 and find(HDIA, "I21")=1 then CHD=1 ;
else if ICD_VERSION=10 and find(HDIA, "I22")=1 then CHD=1;
else if ICD_VERSION=9 and find(HDIA, "410")=1 then CHD=1;
else if ICD_VERSION=9 and find(HDIA, "411B")=1 then CHD=1;
else if ICD_VERSION=8 and find(HDIA, "410")=1 then CHD=1;
else if ICD_VERSION=8 and find(HDIA, "411")=1 then CHD=1;

do i=1 to 3;
if op_a {i} in ("FNG02","FNG05") then CHD=1;
else if  find(op_a {i}, "FNC")=1 then CHD=1;
else if  find(op_a {i}, "FND")=1 then CHD=1;
else if  find(op_a {i}, "FNE")=1 then CHD=1;
else if  find(op_a {i}, "3080")=1 then CHD=1;
else if  find(op_a {i}, "3127")=1 then CHD=1;
else if  find(op_a {i}, "3158")=1 then CHD=1;
end;
if CHD =. then CHD=0;


if ICD_VERSION=10     and find(HDIA, "I63")=1 then IS=1;
else if ICD_VERSION=9 and find(HDIA, "433")=1 then IS=1 ;
else if ICD_VERSION=9 and find(HDIA, "434")=1 then IS=1;
else if ICD_VERSION=8 and find(HDIA, "432")=1 then IS=1;
else if ICD_VERSION=8 and find(HDIA, "433")=1 then IS=1;
else if ICD_VERSION=8 and find(HDIA, "434")=1 then IS=1;
else IS=0;


if CHD=1 or IS=1;

drop OP1-OP3 HDIA  ICD_VERSION i;
run;

* Select events with the checkdate and define preCVD;

proc sort data=cvd; by TWINNR indat; run;
proc sort data=health; by twinnr; run;
data cvd;
merge cvd health (keep=twinnr check);
by twinnr;
if check =. or indat=. then delete;
if indat => check then preCVD=0;
else preCVD=1;
run;

* Subjects with PreCVD (i keep only the last visit befor checkdate);
proc sort data=cvd; by twinnr indat; run;

data preCVD;
set cvd;
by twinnr indat;
if last.twinnr;
keep twinnr preCVD;
where preCVD=1;
run;


* Keep first date after checkdate;
data cvdCHD;
set cvd;
by twinnr indat;
if first.twinnr;
where CHD=1 and preCVD=0;
rename indat=CHD_d;
keep TWINNR OP1-OP3 HDIA indat  CHD;
run;

data cvdIS;
set cvd;
by twinnr indat;
if first.twinnr;
where IS=1 and preCVD=0;
rename indat=IS_d;
keep TWINNR indat IS;
run;

* Merge;
data cvd;
merge cvdCHD cvdIS preCVD;
by twinnr;
run;


*** Cuase of deaths ***;

data cdeath;
set str.link_death_cause;
if MORSAK_NR=1 then output;
run;

*** Deaths ***;

data adeath;
set str.person_info;
deathd=mdy(DEATHMON,DEATHDAY, DEATHYR);
keep TWINNR deathd;
if deathd ne . and deathd <= 18262; * 31DEC2009;
format deathd date9.;
run;


proc sort data=adeath; by twinnr; run;
proc sort data=cdeath; by twinnr; run;

data death;
merge adeath cdeath;
by twinnr;

if find (MORSAK ,"I200") then dCHD=1;
else if find (MORSAK, "411B") then dCHD=1;
else if find(MORSAK, "I21")=1 then dCHD=1 ;
else if find(MORSAK, "I22")=1 then dCHD=1;
else if find(MORSAK, "410")=1 then dCHD=1;
else dCHD=0;


if find(MORSAK, "I63")=1 then dIS=1;
else if find(MORSAK, "433")=1 then dIS=1;
else if find(MORSAK, "434")=1 then dIS=1;
else dIS=0;


* deleting subject with cause of death but not death date;
if deathd =. then delete;
if dCHD ne 1 and dIS ne 1 then dother=1;
else dother=0;

drop MORSAK_NR MORSAK;
run;


* Check no patients checked after death;

proc sort data=death; by TWINNR; run;
proc sort data=health; by twinnr; run;
data death;
merge death health (keep=twinnr check);
by twinnr;
if deathd =. or check=. then delete;
if deathd => check ;
run;


*** Sex and PAIRID***;

data sex;
set twin.V_admin;
keep TWINNR SEX PAIRID;
run;

*** Laboratory data ***;

data labdata;
set twin.V_labdata;
keep TWINNR ANALYSISDESC VALUE;
run;
proc sort data=labdata; by twinnr; run; 
proc transpose data=labdata out=labdata ; by twinnr; var value; id ANALYSISDESC; run;


****** DATA FROM SALT *******;

data smoke;
set salt.X_ibs_az;
keep twinnr XCSMOKER;
run;


****** MERGE *****;
proc sort data=health; by twinnr; run; 
proc sort data=cvd; by twinnr; run; 
proc sort data=death; by twinnr; run; 
proc sort data=sex; by twinnr; run; 
proc sort data=labdata; by twinnr; run; 
proc sort data=smoke; by twinnr; run;

data twingene;
merge health (in=x) cvd death sex labdata smoke;
by twinnr;
if x=1;
run;


*************************************
***** CREATING TWINGENE DATASET *****
*************************************;

data tw.twingene;
set twingene;

* From alphanumeric to numeric;
diabetess=input(DIABETES, best12.);
DIAS_ONE_MMHG_=input(DIAS_ONE_MMHG, best12.);
DIAS_TWO_MMHG_=input(DIAS_TWO_MMHG, best12.);
height=input(HEIGHT_CM, best12.);
SYST_ONE_MMHG_=input(SYST_ONE_MMHG, best12.);
SYST_TWO_MMHG_=input(SYST_TWO_MMHG, best12.);
weight=input(WEIGHT_KG, best12.);
glucos=input(S_Glukos, best12.);
CRP1=input(S_P_CRP_h_gk_nsligt, best12.);
TC=input(S_Kolesterol, best12.);
HDL=input(S_HDL_Kolesterol, best12.);
TG=input(fS_Triglycerider, best12.);
HbA1c=input(B_HbA1c, best12.);
LDL=input(fS_LDL_Kolesterol, best12.);
ApoA1=input(S_Apolipoprotein_A1, best12.);
ApoB=input(S_Apolipoprotein_B, best12.);
Hemoglo=input(B_Hemoglobin, best12.);
CRP2=input(S_P_CRP_h_gk_nsligt2, best12.);

* From alphanumeric to date;
birth_d=mdy(substr(BIRTHDATE,5,2),substr(BIRTHDATE,7,2), substr(BIRTHDATE,1,4));
check_d=check;
death_d=deathd;

* Missing;

if weight in (999,998,997) then weight=.;
if height in (999,998,997) then height =.;
if diabetess in (999,998,997) then diabetess =.;
if preCVD =. then preCVD=0;
if CHD =. then CHD=0;
if IS =. then IS=0;
if dCHD =. then dCHD=0;
if dIS =. then dIS=0;
if dother=. then dother=0;

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


* Outlayer;

if age>110 or age<20 then age=.;
if sbp>250 or sbp<70 then sbp=.;
if dbp>150 or dbp<40 then dbp=.;
if weight>200 or weight<30 then weight=.;
if height>220 or height<100 then height=.;
if height<100 then weight=.;

* Only real checkdate;
if check_d>18262 then delete; * 31DEC2009;


* Setting variables;
if SYST_ONE_MMHG_ ne . and SYST_TWO_MMHG_ ne . then
 sbp=(SYST_ONE_MMHG_+SYST_TWO_MMHG_)/2;
else if SYST_ONE_MMHG_ ne . and SYST_TWO_MMHG_ eq . then 
 sbp=(SYST_ONE_MMHG_);

if DIAS_ONE_MMHG_ ne . and DIAS_ONE_MMHG_ ne . then
 dbp=(DIAS_ONE_MMHG_+DIAS_TWO_MMHG_)/2;
else if DIAS_ONE_MMHG_ ne . and DIAS_ONE_MMHG_ eq . then 
 sbp=(DIAS_ONE_MMHG_);

if CRP1 = . then CRP1=CRP2;
crp=crp1;

smoke=XCSMOKER;

if weight ne . and height ne . then
 bmi=weight/(height/100)**2;
else bmi=.;

age=(check_d-birth_d)/365.25;

* setting outcomes;
if death_d ne . then death=1; else death=0;

*Setting follow-up for survival;
endfup=18262;
  *CHD;
  if death=1 and CHD=0 and dCHD=0 then CHDd=death_d;
  if death=1 and CHD=0 and dCHD=1 then CHDd=death_d;
  if death=1 and CHD=1 and dCHD=1 then CHDd=CHD_d;
  if death=1 and CHD=1 and dCHD=0 then CHDd=CHD_d;
  if death=0 and CHD=1 and dCHD=0 then CHDd=CHD_d;
  if death=0 and CHD=0 and dCHD=0 then  CHDd=endfup;
  if CHD =1 or dCHD=1 then CHD=1;
  survCHD=CHDd-check_d;
 *IS;
  if death=1 and IS=0 and dIS=0 then ISd=death_d;
  if death=1 and IS=0 and dIS=1 then ISd=death_d;
  if death=1 and IS=1 and dIS=1 then ISd=IS_d;
  if death=1 and IS=1 and dIS=0 then ISd=IS_d;
  if death=0 and IS=1 and dIS=0 then ISd=IS_d;
  if death=0 and IS=0 and dIS=0 then  ISd=endfup;
  if IS =1 or dIS=1 then IS=1;
  survIS=ISd-check_d;


* Label;
label

age='Age'
sex='Sex'
diabetess='Anamnestic diabetes'
height='Height in cm.'
weight='Weight in Kg.'
betab='Beta blockers'
ACE_angio='ACE/AngiotensinII inibithors'
calcb='Calcium blockers'
diur='Diuretics'
statins='Statins'
fibr='Fibrates'
preCVD = 'Previous CVD'
CHD = 'CHD'
IS = 'Ischemic Stroke'
glucos='Glucose mmol/L'
crp='CRP mg/L'
TC='Total Cholesterol mmol/L'
HDL='HDL Cholesterol mmol/L'
TG='Triglyceride mmol/L'
HbA1c='HbA1c %'
LDL='LDL  Cholesterol mmol/L'
ApoA1='Apolipoprotein A-I g/L'
ApoB='Apolipoprotein B g/L'
Hemoglo='Hemoglobine g/L'
sbp='Systolic Blood Pressure mm/Hg'
dbp='Diastolic Blood Pressure mm/Hg'
birth_d='Birth date'
check_d='Check date'
death_d='Death date (all causes)'
smoke='Current smoker'
bmi='Body Max Index'
death='Death/alive';


drop DIABETES DIAS_ONE_MMHG DIAS_TWO_MMHG  HEIGHT_CM 
     SYST_ONE_MMHG SYST_TWO_MMHG WEIGHT_KG S_Glukos S_P_CRP_h_gk_nsligt S_Kolesterol
	 S_HDL_Kolesterol fS_Triglycerider B_HbA1c fS_LDL_Kolesterol S_Apolipoprotein_A1 S_Apolipoprotein_B
     B_Hemoglobin S_P_CRP_h_gk_nsligt2 BIRTHDATE deathd  _NAME_ _LABEL_ check
     DIAS_ONE_MMHG_ DIAS_TWO_MMHG_ SYST_ONE_MMHG_ SYST_TWO_MMHG_ CRP1 CRP2 XCSMOKER 
     dCHD dIS CHD_d IS_d dother endfup HART_MEDICINE_TYPE;

format birth_d date9. check_d date9. death_d date9. endfup date9. CHDd date9. ISd date9.; 
run;


*********************************************** 
***** TWINGENE WITH UNRELATED INDIVIDUALS *****
***********************************************;

proc sort data=tw.twingene; by PAIRID; run;

proc surveyselect data=tw.twingene sampsize=1 out=tw.twingene_unr seed=333 method=srs;
strata PAIRID; 
run;


******************************
***** PREPARE AND EXPORT *****
******************************;

* Create some important variables;

data wc;
set tw.twingene_unr;

*we put together CHD and IS and we call it CVD;
if CHD=1 or IS=1 then CVD=1;
else CVD=0;
if survCHD <= survIS then survCVD=survCHD;
else survCVD=survIS;

*we want the CVD date;

if CHDd <= ISd then CVDd=CHDd;
else CVDd=ISd;

* Start and stop, as everyone enter at time 0;
start=0;
stop=survCVD;

* Define anti-hypertensive treatmenr;

if betab=1 or ACE_angio=1 or calcb=1 or diur=1 then antihyp=1;
else antihyp=0;

*Merging indicator;
a=1;

*Delete previous CVD (or CHD or IS);
if preCVD=1 then delete;

* Delete missing data;

if sex=. or TC=. or HDL=. or age=. or sbp=. or smoke=. or antihyp=. or apoA1= . or diabetess= . then delete;

*Keeping usefull varuiables;

keep twinnr sex TC HDL age sbp smoke antihyp CVD start stop a check_d CVDd apoA1 diabetess;

run;


* Breaking ties;
* REMEMBER THAT THIS SYSTEM WORKS ONLY WITH MAX 2 EQUAL FAILURE TIMES,
  IF THREE OR MORE DO IT MANUALLY!!!;

proc freq data=wc noprint; table stop / out=freq; where cvd=1; run; 
proc sort data=freq; by stop; run;
proc sort data=wc; by stop descending CVD; run;
data wc; merge wc freq; by stop; 
 if COUNT > 1 and cvd=1 then do;
    if first.stop then stop=stop+0.1;
 end;
 if COUNT = 3 then do;
  if first.stop then stop=stop+0.1;
  else if last.stop then stop=stop-0.1;
 end;
drop PERCENT COUNT;
run;
proc sort data=wc; by stop; run;

* Crp is standardized;
*proc standard data=wc mean=0 std=1 out=wc; var crp; run;

*Delete temporary datasets;
proc datasets library = work nolist;
delete freq twingene health cvd death sex labdata smoke;
quit;

**** EXPORT FOR R PROGRAMS;

PROC EXPORT DATA= WORK.Wc 
            OUTFILE= "&path.\Data\wc.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




