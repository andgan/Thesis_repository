/*----------------------------------------------
Filename: import_GOSH.sas
Study: genscore
Author: Andrea Ganna
Date: 25JUN2011
Updated: 20SEP2011 - Added lipid lowering drougs (no one has it!)
         20OCT2011 - Keep only SATSA-3
         02DEC2011 - Added gender IPT1
		 03DEC2011 - Update IPT3 for satsa		
         05DEC2011 - Now SATSA 3 have complete info about diabetes and harmony about smoke
Purpose: Import data from 4 different study (GOSH): Gender Octo Satsa Harmony tobe use in the genscore project
Note:
-----------------------------------------------
Data used: databases from ORACLE and from P:/
Data created: gosh.GOSH
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/


**** ASSIGN LIBRARIES;

*%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;
%let path=E:\Phd_KI\Genscore;


libname gosh "&path.\Data\GOSH\Final";


***************************************
************  GENDER  *****************
***************************************;
** APO A-1 used instead of HDL-C

* Blood and variables from IPT1, variables from IPT2;
/*
data age_sex;
set gender.INCIDENT_COMPLETE_080812;
keep TWINNR PAIRID check_g;
check_g=mdy(substr(G_date1 ,6,2),substr(G_date1 ,9,2), substr(G_date1 ,1,4));
format check_g date9. ;
run;

data lab;
set gender.GENDANALYS_LAB;
keep TWINNR CHOL APOA1;
run;


data FHS;
set gender.G2A;
sbp_g=mean(VA220,VA222);
weight=VA230;
height=VA231;
bmi_g=weight/(height/100)**2;
diab_g=VA273;

if substr(UPCASE(VA300),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA375),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA380),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA385),1,3)in ("C02","C03","C07","C08","C09") 
   or  substr(UPCASE(VA390),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA395),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA400),1,3)in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(VA405),1,3)in ("C02","C03","C07","C08","C09")
   or  substr(UPCASE(VA410),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA415),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(VA420),1,3)in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(VA425),1,3)in ("C02","C03","C07","C08","C09") then antihyp_g =1;
   else antihyp_g =0;

if substr(UPCASE(VA300),1,3)in ("C10") or  substr(UPCASE(VA375),1,3)="C10" or  substr(UPCASE(VA380),1,3)in ("C10") or  substr(UPCASE(VA385),1,3)in ("C10") 
   or  substr(UPCASE(VA390),1,3)in ("C10") or  substr(UPCASE(VA395),1,3)in ("C10") or  substr(UPCASE(VA400),1,3)in ("C10")  or  substr(UPCASE(VA405),1,3)in ("C10")
   or  substr(UPCASE(VA410),1,3)in ("C10") or  substr(UPCASE(VA415),1,3)in ("C10") or  substr(UPCASE(VA420),1,3)in ("C10")  or  substr(UPCASE(VA425),1,3)in ("C10") then drug_lipid_g =1;
   else drug_lipid_g =0;

if substr(UPCASE(VA300),1,3)="A10" or  substr(UPCASE(VA375),1,3)="A10" or  substr(UPCASE(VA380),1,3)="A10" or  substr(UPCASE(VA385),1,3)="A10" 
   or  substr(UPCASE(VA390),1,3)="A10" or  substr(UPCASE(VA395),1,3)="A10" or  substr(UPCASE(VA400),1,3)="A10"  or  substr(UPCASE(VA405),1,3)="A10"
   or  substr(UPCASE(VA410),1,3)="A10" or  substr(UPCASE(VA415),1,3)="A10" or  substr(UPCASE(VA420),1,3)="A10"  or  substr(UPCASE(VA425),1,3)="A10" then diab_treat =1;
   else diab_treat =0;


keep TWINNR sbp_g diab_g bmi_g antihyp_g diab_treat drug_lipid_g;
run;

data smoke;
set gender.G2B;
if VA450 = 1 or VA450=3 then smoke_g=1;
else smoke_g=0;
keep TWINNR smoke_g;
run;


proc sort data=age_sex; by TWINNR; run;
proc sort data=lab; by TWINNR; run;
proc sort data=FHS; by TWINNR; run;
proc sort data=smoke; by TWINNR; run;


data gender;
merge age_sex  FHS  lab smoke ;
by TWINNR;

*Indicator;
gender=1;

* character to numeric;
tc_g=input(CHOL,best12.);
apoa_g=input(APOA1,best12.);

* Calculate diabetes;
if diab_g = . and diab_treat = . then diab_g = .;
else if diab_g =1 or diab_treat=1 then diab_g = 1;
else diab_g = 0;

* Keep only people wih first visit;
if check_g ne .;

* MISSING;
sbp_m=missing(sbp_g);
bmi_m=missing(bmi_g);
diab_m=missing(diab_g);
antihyp_m=missing(antihyp_g);
drug_lipid_m=missing(drug_lipid_g);
smoke_m=missing(smoke_g);
tc_m=missing(tc_g);
apoa_m=missing(apoa_g);

miss_g=sbp_m+bmi_m+diab_m+antihyp_m+drug_lipid_m+smoke_m+tc_m+apoa_m+drug_lipid_m;


drop diab_treat CHOL APOA1 sbp_m bmi_m diab_m antihyp_m drug_lipid_m smoke_m tc_m apoa_m;
run;*/



*** VARIABLES FROM IPT1 ****;
libname gipt1 "&path.\Data\GOSH";

data age_sex;
set gender.INCIDENT_COMPLETE_080812;
keep TWINNR PAIRID check_g;
check_g=mdy(substr(G_date1 ,6,2),substr(G_date1 ,9,2), substr(G_date1 ,1,4));
format check_g date9. ;
run;

data FHS;
set Gipt1.Genderorginal_andrea;
sbp_g=mean(VA220,VA222);
weight=VA230;
height=VA231;
bmi_g=weight/(height/100)**2;
diab_g=VA273;

if VA450 = 1 or VA450=3 then smoke_g=1;
else smoke_g=0;

if substr(UPCASE(lb1),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb2),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb3),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb4),1,3)in ("C02","C03","C07","C08","C09") 
   or  substr(UPCASE(lb5),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb6),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb6),1,3)in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(lb8),1,3)in ("C02","C03","C07","C08","C09")
   or  substr(UPCASE(lb9),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb10),1,3)in ("C02","C03","C07","C08","C09") or  substr(UPCASE(lb11),1,3)in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(lb12),1,3)in ("C02","C03","C07","C08","C09") then antihyp_g =1;
   else antihyp_g =0;

if substr(UPCASE(lb1),1,3)in ("C10") or  substr(UPCASE(lb2),1,3)="C10" or  substr(UPCASE(lb3),1,3)in ("C10") or  substr(UPCASE(lb4),1,3)in ("C10") 
   or  substr(UPCASE(lb5),1,3)in ("C10") or  substr(UPCASE(lb6),1,3)in ("C10") or  substr(UPCASE(lb7),1,3)in ("C10")  or  substr(UPCASE(lb8),1,3)in ("C10")
   or  substr(UPCASE(lb9),1,3)in ("C10") or  substr(UPCASE(lb10),1,3)in ("C10") or  substr(UPCASE(lb11),1,3)in ("C10")  or  substr(UPCASE(lb12),1,3)in ("C10") then drug_lipid_g =1;
   else drug_lipid_g =0;

if substr(UPCASE(lb1),1,3)="A10" or  substr(UPCASE(lb2),1,3)="A10" or  substr(UPCASE(lb3),1,3)="A10" or  substr(UPCASE(lb4),1,3)="A10" 
   or  substr(UPCASE(lb5),1,3)="A10" or  substr(UPCASE(lb6),1,3)="A10" or  substr(UPCASE(lb7),1,3)="A10"  or  substr(UPCASE(lb8),1,3)="A10"
   or  substr(UPCASE(lb9),1,3)="A10" or  substr(UPCASE(lb10),1,3)="A10" or  substr(UPCASE(lb11),1,3)="A10"  or  substr(UPCASE(lb12),1,3)="A10" then diab_treat =1;
   else diab_treat =0;


tc_g=kol;
apoa_g=apoa;

keep TWINNR sbp_g diab_g bmi_g antihyp_g diab_treat drug_lipid_g tc_g apoa_g smoke_g;
run;



proc sort data=age_sex; by TWINNR; run;
proc sort data=FHS; by TWINNR; run;

data gender;
merge age_sex FHS;
by TWINNR;

*Indicator;
gender=1;

* Calculate diabetes;
if diab_g = . and diab_treat = . then diab_g = .;
else if diab_g =1 or diab_treat=1 then diab_g = 1;
else diab_g = 0;

* Keep only people wih first visit;
if check_g ne .;

* MISSING;
sbp_m=missing(sbp_g);
bmi_m=missing(bmi_g);
diab_m=missing(diab_g);
antihyp_m=missing(antihyp_g);
drug_lipid_m=missing(drug_lipid_g);
smoke_m=missing(smoke_g);
tc_m=missing(tc_g);
apoa_m=missing(apoa_g);

miss_g=sbp_m+bmi_m+diab_m+antihyp_m+drug_lipid_m+smoke_m+tc_m+apoa_m+drug_lipid_m;


drop diab_treat CHOL APOA1 sbp_m bmi_m diab_m antihyp_m drug_lipid_m smoke_m tc_m apoa_m;
run;



***************************************
************  OCTO  *******************
***************************************;

* blood has been withdrawn at IPT 2 (wave 2), but we use risk factors from IPT1;

data age_sex;
set octo.OCTO_INCIDENT_080930;
keep TWINNR PAIRID check_o;
check_o=mdy(substr(O_date1,6,2),substr(O_date1,9,2), substr(O_date1,1,4));
format check_o date9. ;
run;

data lab;
set octo.BLOOD_RESULTS;
keep TWINNR KOLX HDLX;
run;

data smoke;
set octo.WAVE1;
if V473=. then smoke_o=.;
else if V473 = 1 or V473=3 then smoke_o =1 ;
else smoke_o =0;
keep TWINNR smoke_o;
run;

data FHS;
set octo.WAVE1;
sbp_o=V232;
weight=V244;
height=V246;
diab_o=V284;
bmi_o=weight/(height/100)**2;
keep TWINNR sbp_o diab_o bmi_o  ;
run;



data antihyp_diab_t;
set octo.WAVE1;
if substr(UPCASE(V401),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V406),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V411),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V416),1,3) in ("C02","C03","C07","C08","C09") 
   or  substr(UPCASE(V421),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V426),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V431),1,3) in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(V436),1,3) in ("C02","C03","C07","C08","C09")
   or  substr(UPCASE(V441),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V446),1,3) in ("C02","C03","C07","C08","C09") or  substr(UPCASE(V451),1,3) in ("C02","C03","C07","C08","C09")  or  substr(UPCASE(V456),1,3) in ("C02","C03","C07","C08","C09") then antihyp_o =1;
   else antihyp_o =0;

if substr(UPCASE(V401),1,3)="C10" or  substr(UPCASE(V406),1,3)="C10" or  substr(UPCASE(V411),1,3)="C10" or  substr(UPCASE(V416),1,3)="C10" 
   or  substr(UPCASE(V421),1,3)="C10" or  substr(UPCASE(V426),1,3)="C10" or  substr(UPCASE(V431),1,3)="C10"  or  substr(UPCASE(V436),1,3)="C10"
   or  substr(UPCASE(V441),1,3)="C10" or  substr(UPCASE(V446),1,3)="C10" or  substr(UPCASE(V451),1,3)="C10"  or  substr(UPCASE(V456),1,3)="C10" then drug_lipid_o =1;
   else drug_lipid_o =0;

if substr(UPCASE(V401),1,3)="A10" or  substr(UPCASE(V406),1,3)="A10" or  substr(UPCASE(V411),1,3)="A10" or  substr(UPCASE(V416),1,3)="A10" 
   or  substr(UPCASE(V421),1,3)="A10" or  substr(UPCASE(V426),1,3)="A10" or  substr(UPCASE(V431),1,3)="A10"  or  substr(UPCASE(V436),1,3)="A10"
   or  substr(UPCASE(V441),1,3)="A10" or  substr(UPCASE(V446),1,3)="A10" or  substr(UPCASE(V451),1,3)="A10"  or  substr(UPCASE(V456),1,3)="A10" then diab_treat =1;
   else diab_treat =0;

keep TWINNR antihyp_o diab_treat drug_lipid_o;
run;



*** MERGE ***;


proc sort data=age_sex; by TWINNR; run;
proc sort data=lab; by TWINNR; run;
proc sort data=FHS; by TWINNR; run;
proc sort data=smoke; by TWINNR; run;
proc sort data=antihyp_diab_t; by TWINNR; run;



data octo;
merge age_sex lab  smoke FHS  antihyp_diab_t ;
by TWINNR;

*Indicator;
octo=1;

* Calculate diabetes;
if diab_o = . and diab_treat ne 1 then diab_o = .;
else if diab_o =1 or diab_treat=1 then diab_o = 1;
else diab_o = 0;

* keep only subjects in this second visit;
if check_o ne . ;

* Rename;
tc_o=KOLX;
hdl_o=HDLX;

* MISSING;
sbp_m=missing(sbp_o);
bmi_m=missing(bmi_o);
diab_m=missing(diab_o);
antihyp_m=missing(antihyp_o);
drug_lipid_m=missing(drug_lipid_o);
smoke_m=missing(smoke_o);
tc_m=missing(tc_o);
hdl_m=missing(hdl_o);

miss_o=sbp_m+bmi_m+diab_m+antihyp_m+smoke_m+tc_m+hdl_m+drug_lipid_m;

* Keep only people with visit;
if check_o ne .;

drop  diab_treat KOLX HDLX sbp_m bmi_m diab_m antihyp_m smoke_m tc_m hdl_m drug_lipid_m;
run;



*****************************************
************  SATSA_1  ******************
*****************************************;



data age_sex_1;
set satsa.SATSA_INCIDENT_110207;
keep TWINNR check_s_1 ;
check_s_1=mdy(substr(S_DATE1 ,6,2),substr(S_DATE1 ,9,2), substr(S_DATE1 ,1,4));
format check_s_1 date9. ;
run;

* set libaname to get data directely from P: because not in oracle;
libname temp1 v6 "P:\twins\str_projects\internal_base\SATSA\ipt1";

data lab_1;
set temp1.Ipt1bld2;
keep IDA REGA CHOLEST HDL;
run;


data BMI_sbp_1;
set satsa.SATSA_IPT1VER3;
sbp_s_1=BPSRES;
weight= WEIGHT;
height= LENGTH;
bmi_s_1=weight/(height/100)**2;
keep TWINNR sbp_s_1  bmi_s_1;
run;


data smoke_1;
set satsa.SATSA_Q1_REDSCALE;

if SMOKSTA = . then smoke_s_1 = .;
else if SMOKSTA = 3 then smoke_s_1 = 1;
else smoke_s_1 = 0; 

keep TWINNR smoke_s_1;
run;


* set libaname to get data directely from P: because not in oracle;
libname temp2  "P:\twins\str_projects\internal_base\SATSA\ipt1";
data treat_1;
set temp2.Depharm;
if CLASSI=1 or CLASSJ=1 or  BETA=1 or CALANT=1 or ACE=1 then antihyp_s_1=1;
if CLASSO=1 then diab_treat=1;
*keep IDA REGA antihyp_s_1 diab_treat;
run;


data diab1_1;
set satsa.SATSA_Q1_RED;

if DIABETE = . then diab1 = .;
else if DIABETE = 2 then diab1 = 1;
else diab1 = 0; 

keep TWINNR IDA  REGA diab1 ;
run;

data diab2_1;
set temp2.Ipt1ill1;
diab2=MDIABETE;
keep IDA REGA diab2;
run;


* Previous merge to get rid of IDA;
proc sort data=treat_1; by REGA IDA; run;
proc sort data=diab1_1; by REGA IDA; run;
proc sort data=diab2_1; by REGA IDA; run;
proc sort data=lab_1; by REGA IDA; run;


data mm_1;
merge treat_1 diab1_1 diab2_1 lab_1;
by REGA IDA;

* Calculate diabetes;
if diab1 = . and diab_treat = . then diab_s_1 = .;
else if diab1 =1 or diab2=1 or diab_treat=1 then diab_s_1 = 1;
else diab_s_1 = 0;

if antihyp_s_1=. then antihyp_s_1=0;
drug_lipid_s_1=0;

* Some subjects are without TWINNR;
if TWINNR ne .;

keep TWINNR antihyp_s_1 diab_s_1 CHOLEST HDL drug_lipid_s_1;

run;



*** MERGE ***;

proc sort data=age_sex_1; by TWINNR; run;
proc sort data=mm_1; by TWINNR; run;
proc sort data=BMI_sbp_1; by TWINNR; run;
proc sort data=smoke_1; by TWINNR; run;

** satsa from IPT1;

data satsa_1;
merge age_sex_1 (in=x) mm_1  BMI_sbp_1 smoke_1 ;
by TWINNR;

*Indicator;
satsa_1=1;

* character to numeric;
tc_s_1= CHOLEST;
hdl_s_1= HDL;
* Keep only people with visit;
if check_s_1 ne .;

* MISSING;
sbp_m=missing(sbp_s_1);
bmi_m=missing(bmi_s_1);
diab_m=missing(diab_s_1);
antihyp_m=missing(antihyp_s_1);
drug_lipid_m=missing(drug_lipid_s_1);
smoke_m=missing(smoke_s_1);
tc_m=missing(tc_s_1);
hdl_m=missing(hdl_s_1);

miss_s_1=sbp_m+bmi_m+diab_m+antihyp_m+smoke_m+tc_m+hdl_m+drug_lipid_m;


drop  CHOLEST HDL sbp_m bmi_m diab_m antihyp_m smoke_m tc_m hdl_m drug_lipid_m;
run;



*****************************************
************  SATSA_2  ******************
*****************************************;



data age_sex_2;
set satsa.SATSA_INCIDENT_110207;
keep TWINNR check_s_2 ;
check_s_2=mdy(substr(S_DATE2 ,6,2),substr(S_DATE2 ,9,2), substr(S_DATE2 ,1,4));
format check_s_2 date9. ;
run;

* set libaname to get data directely from P: because not in oracle;
libname temp3 "P:\twins\str_projects\internal_base\SATSA\ipt2";

data lab_2;
set temp3.Ipt2bld2;
keep IDA REGA QCHOLEST QAPOA;
run;


data BMI_sbp_2;
set satsa.SATSA_IPT2E;
sbp_s_2=QBPSRES;
weight= QWEIGHT;
height= QLENGTH;
bmi_s_2=weight/(height/100)**2;
keep TWINNR sbp_s_2  bmi_s_2;
run;


data treat_2;
set temp3.ipt2med (rename=(ida=IDA_ rega=REGA_));
if substr(UPCASE(FASS98),1,3)in ("C02","C03","C07","C08","C09") then antihyp_s_2=1;
if substr(UPCASE(FASS98),1,3)in ("C10") then drug_lipid_s_2=1;
if substr(UPCASE(FASS98),1,3)="A10" then diab_treat=1;
if antihyp_s_2=1 or diab_treat=1 or drug_lipid_s_2=1;

* This is to give a best12. format comaptible with the other datasets;
REGA=input(REGA_,best12.);
IDA=input(IDA_,best12.);
keep IDA REGA antihyp_s_2 diab_treat drug_lipid_s_2;
run;

proc sort data=treat_2; by REGA IDA; run;
data diab_treat_2;
set treat_2;
by REGA IDA;
if first.IDA;
where diab_treat ne .;
drop antihyp_s_2 drug_lipid_s_2;
run;

data antihyp_2;
set treat_2;
by REGA IDA;
if first.IDA;
where antihyp_s_2 ne .;
drop diab_treat drug_lipid_s_2;
run;

data drug_lipid_2;
set treat_2;
by REGA IDA;
if first.IDA;
where drug_lipid_s_2 ne .;
drop diab_treat antihyp_s_2;
run;


* set libaname to get data directely from P: because not in oracle;
libname temp4 v604 "P:\twins\str_projects\internal_base\SATSA\ipt2";

data smoke_2;
set temp4.Qscale1;
if QSMOKSTA=3 then smoke_s_2=1; else smoke_s_2=0;
keep IDA REGA smoke_s_2;
run;

data diab_2;
set satsa.SATSA_IPT2Q;

if QDIABETE = . then diab = .;
else if QDIABETE = 2 then diab = 1;
else diab = 0; 

keep TWINNR IDA REGA diab ;
run;


* Previous merge to get rid of IDA;
proc sort data=lab_2; by REGA IDA; run;
proc sort data=diab_treat_2; by REGA IDA; run;
proc sort data=antihyp_2; by REGA IDA; run;
proc sort data=drug_lipid_2; by REGA IDA; run;
proc sort data=diab_2; by REGA IDA; run;
proc sort data=smoke_2; by REGA IDA; run;


data mm_2;
merge lab_2 diab_treat_2 antihyp_2 drug_lipid_2 smoke_2 diab_2;
by REGA IDA;

* Calculate diabetes;
if diab = . and diab_treat = . then diab_s_2 = .;
else if diab =1 or diab_treat=1 then diab_s_2 = 1;
else diab_s_2 = 0;

if antihyp_s_2=. then antihyp_s_2=0;
if drug_lipid_s_2=. then drug_lipid_s_2=0;
* Some subjects are without TWINNR;
if TWINNR ne .;

keep TWINNR antihyp_s_2 diab_s_2 smoke_s_2 QCHOLEST QAPOA drug_lipid_s_2;
run;



*** MERGE ***;

proc sort data=age_sex_2; by TWINNR; run;
proc sort data=mm_2; by TWINNR; run;
proc sort data=BMI_sbp_2; by TWINNR; run;


** satsa from IPT2;

data satsa_2;
merge age_sex_2 (in=x) mm_2  BMI_sbp_2 ;
by TWINNR;

*Indicator;
satsa_2=1;

* character to numeric;
tc_s_2= QCHOLEST;
hdl_s_2= QAPOA;
* Keep only people with visit;
if check_s_2 ne .;

* MISSING;
sbp_m=missing(sbp_s_2);
bmi_m=missing(bmi_s_2);
diab_m=missing(diab_s_2);
antihyp_m=missing(antihyp_s_2);
drug_lipid_m=missing(drug_lipid_s_2);
smoke_m=missing(smoke_s_2);
tc_m=missing(tc_s_2);
hdl_m=missing(hdl_s_2);

miss_s_2=sbp_m+bmi_m+diab_m+antihyp_m+smoke_m+tc_m+hdl_m+drug_lipid_m;


drop  QCHOLEST QAPOA sbp_m bmi_m diab_m antihyp_m smoke_m tc_m hdl_m drug_lipid_m;
run;



*****************************************
************  SATSA_3  ******************
*****************************************;

** DNA from SATSA_3;
** APOA1 is used instead of HDL-C;

data age_sex_3;
set satsa.SATSA_INCIDENT_110207;
keep TWINNR check_s_3;
check_s_3=mdy(substr(S_DATE3,6,2),substr(S_DATE3,9,2), substr(S_DATE3,1,4));
format check_s_3 date9. ;
run;


data BMI_sbp_3;
set satsa.SATSA_IPT3;
sbp_s_3=XBPCRES;
weight=XWEIGHT;
height=XLENGTH;
bmi_s_3=weight/(height/100)**2;
keep TWINNR sbp_s_3  bmi_s_3;
run;


* set libaname to get data directely from P: because not in oracle;
libname temp5 v6 "P:\twins\str_projects\internal_base\SATSA\ipt3";

data lab_3;
set temp5.ipt3blod;
keep IDA REGA CHOL APOA1 ;
run;

data treat_3;
set temp5.ipt3med;
if substr(UPCASE(FASS98),1,3)in ("C02","C03","C07","C08","C09") then antihyp_s_3=1;
if substr(UPCASE(FASS98),1,3)in ("C10") then drug_lipid_s_3=1;
if substr(UPCASE(FASS98),1,3)="A10" then diab_treat=1;
if antihyp_s_3=1 or diab_treat=1;
keep IDA REGA antihyp_s_3 diab_treat drug_lipid_s_3;
run;

proc sort data=treat_3; by REGA IDA; run;
data diab_treat_3;
set treat_3;
by REGA IDA;
if first.IDA;
where diab_treat ne .;
drop antihyp_s_3 drug_lipid_s_3;
run;

data antihyp_3;
set treat_3;
by REGA IDA;
if first.IDA;
where antihyp_s_3 ne .;
drop diab_treat drug_lipid_s_3;
run;

data drug_lipid_3;
set treat_3;
by REGA IDA;
if first.IDA;
where drug_lipid_s_3 ne .;
drop diab_treat antihyp_s_3;
run;

* set libaname to get data directely from P: because not in oracle;
libname temp6 v604 "P:\twins\str_projects\internal_base\SATSA\ipt3";

data smoke_3;
set temp6.Xscale;
if XSMOKSTA=3 then smoke_s_3=1; else smoke_s_3=0;
keep REGA IDA smoke_s_3;
run;

** DIABETES;

data diab_3a;
set satsa.SATSA_IPT3Q;

if XDIABETE = . then diab1 = .;
else if XDIABETE = 2 then diab1 = 1;
else diab1 = 0; 

keep diab1 TWINNR ;
run;


* Some subjects have missing diabetes from IPT3Q, get it from Q4;
data diab_3b;
set satsa.SATSA_Q4;

if HDIABETE = . then diab2 = .;
else if HDIABETE = 2 then diab2 = 1;
else diab2 = 0; 

keep TWINNR diab2;
run;

data admino;
set satsa.SATSA_ADMIN_RESPONSE_0406;
if IPT3=1;
keep TWINNR REGA IDA IPT3;
run;


* First we merge diabetes with admino to get the right diabetes;
proc sort data=diab_3a; by TWINNR; run;
proc sort data=diab_3b; by TWINNR; run;
proc sort data=admino; by TWINNR; run;

data adminodiab;
merge admino diab_3a diab_3b;
by TWINNR;
diab=diab1;
if diab1 = . then diab=diab2;
keep TWINNR REGA IDA IPT3 diab;
run;


* Previous merge to get rid of IDA;
proc sort data=lab_3; by REGA IDA; run;
proc sort data=diab_treat_3; by REGA IDA; run;
proc sort data=antihyp_3; by REGA IDA; run;
proc sort data=drug_lipid_3; by REGA IDA; run;
proc sort data=smoke_3; by REGA IDA; run;
proc sort data=adminodiab; by REGA IDA; run;


data mm_3;
merge lab_3 diab_treat_3 antihyp_3 drug_lipid_3 smoke_3 adminodiab;
by REGA IDA;

* Calculate diabetes;
if diab = . and diab_treat = . then diab_s_3 = .;
else if diab =1 or diab_treat=1 then diab_s_3 = 1;
else diab_s_3 = 0;

if antihyp_s_3=. then antihyp_s_3=0;
if drug_lipid_s_3=. then drug_lipid_s_3=0;
* Some subjects are without TWINNR;
if TWINNR ne .;

keep TWINNR antihyp_s_3 drug_lipid_s_3 diab_s_3 smoke_s_3 CHOL APOA1;
run;



*** MERGE ***;

proc sort data=age_sex_3; by TWINNR; run;
proc sort data=mm_3; by TWINNR; run;
proc sort data=BMI_sbp_3; by TWINNR; run;


** satsa from IPT3;

data satsa_3;
merge age_sex_3 (in=x) mm_3  BMI_sbp_3 ;
by TWINNR;

*Indicator;
satsa_3=1;

* character to numeric;
tc_s_3=CHOL;
hdl_s_3=APOA1;
* Keep only people with visit;
if check_s_3 ne .;

* MISSING;
sbp_m=missing(sbp_s_3);
bmi_m=missing(bmi_s_3);
diab_m=missing(diab_s_3);
antihyp_m=missing(antihyp_s_3);
drug_lipid_m=missing(drug_lipid_s_3);
smoke_m=missing(smoke_s_3);
tc_m=missing(tc_s_3);
hdl_m=missing(hdl_s_3);

miss_s_3=sbp_m+bmi_m+diab_m+antihyp_m+smoke_m+tc_m+hdl_m+drug_lipid_m;


drop CHOL APOA1 sbp_m bmi_m diab_m antihyp_m smoke_m tc_m hdl_m drug_lipid_m;
run;


*****************************************
************  HARMONY  ******************
*****************************************;


data FHS;
set harmony.HARMONY_SOMATIC;
keep TWINNR check_h sbp_h bmi_h;
if BLOOD_PRES_REST_PRESSURE_V = "SYSTOL 160, MANUELLT" then BLOOD_PRES_REST_PRESSURE_V = "160";
sbp_h=input(scan(BLOOD_PRES_REST_PRESSURE_V,1,'/'), best12.);
weight=HEIGHT_WEIGHT_WEIGHT_VALUE;
height=HEIGHT_WEIGHT_HEIGHT_VALUE;
bmi_h=weight/(height/100)**2;
check_h=datepart(DATUM_VALUE);
format check_h date9.;
run;


data lab;
set harmony.TWISST_BLOODRESULT;
keep TWINNR KOL HDL HBA1C;

*Delete some suplicated subjects;
if TWINNR=5601 and FASTING_HOURS="9" then delete;
if TWINNR=140192 and FASTING_HOURS="8" then delete;
if TWINNR=171091 and FASTING_HOURS="9" then delete;

where project="Harmony";
run;


* set libaname to get data directely from P: because on oracle is not working;
libname temp2 "P:\twins\str_projects\internal_base\Harmony\data";
data diab;
set temp2.harmony_medhx1;

if diabetes = 998 or diabetes=3 then diab_h=.;
else if diabetes = 0 then diab_h=1;
else diab_h=0;

keep TWINNR diab_h;
run;

data diab_treat_antihyp;
set temp2.harmony_medhx2;

array antihyp_d {20} 
lakemedel_1_funktion_Antihyperte 
lakemedel_2_funktion_Antihyperte 
lakemedel_3_funktion_Antihyperte
lakemedel_4_funktion_Antihyperte
lakemedel_5_funktion_Antihyperte
lakemedel_6_funktion_Antihyperte
lakemedel_7_funktion_Antihyperte
lakemedel_8_funktion_Antihyperte
lakemedel_9_funktion_Antihyperte
lakemedel_10_funktion_Antihypert
lakemedel_11_funktion_Antihypert
lakemedel_12_funktion_Antihypert
lakemedel_13_funktion_Antihypert
lakemedel_14_funktion_Antihypert
lakemedel_15_funktion_Antihypert
lakemedel_16_funktion_Antihypert
lakemedel_17_funktion_Antihypert
lakemedel_18_funktion_Antihypert
lakemedel_19_funktion_Antihypert
lakemedel_20_funktion_Antihypert;

array diuretic_d {20} 
lakemedel_1_funktion_Diretika
lakemedel_2_funktion_Diretika
lakemedel_3_funktion_Diretika
lakemedel_4_funktion_Diretika
lakemedel_5_funktion_Diretika
lakemedel_6_funktion_Diretika
lakemedel_7_funktion_Diretika
lakemedel_8_funktion_Diretika
lakemedel_9_funktion_Diretika
lakemedel_10_funktion_Diretika
lakemedel_11_funktion_Diretika
lakemedel_12_funktion_Diretika
lakemedel_13_funktion_Diretika
lakemedel_14_funktion_Diretika
lakemedel_15_funktion_Diretika
lakemedel_16_funktion_Diretika
lakemedel_17_funktion_Diretika
lakemedel_18_funktion_Diretika
lakemedel_19_funktion_Diretika
lakemedel_20_funktion_Diretika;

array insulin_d {20} 
lakemedel_1_funktion_Insulin 
lakemedel_2_funktion_Insulin 
lakemedel_3_funktion_Insulin
lakemedel_4_funktion_Insulin
lakemedel_5_funktion_Insulin
lakemedel_6_funktion_Insulin
lakemedel_7_funktion_Insulin
lakemedel_8_funktion_Insulin
lakemedel_9_funktion_Insulin
lakemedel_10_funktion_Insulin
lakemedel_11_funktion_Insulin
lakemedel_12_funktion_Insulin
lakemedel_13_funktion_Insulin
lakemedel_14_funktion_Insulin
lakemedel_15_funktion_Insulin
lakemedel_16_funktion_Insulin
lakemedel_17_funktion_Insulin
lakemedel_18_funktion_Insulin
lakemedel_19_funktion_Insulin
lakemedel_20_funktion_Insulin;

array perorala_d {20} 
lakemedel_1_funktion_Perorala_an
lakemedel_2_funktion_Perorala_an
lakemedel_3_funktion_Perorala_an
lakemedel_4_funktion_Perorala_an
lakemedel_5_funktion_Perorala_an
lakemedel_6_funktion_Perorala_an
lakemedel_7_funktion_Perorala_an
lakemedel_8_funktion_Perorala_an
lakemedel_9_funktion_Perorala_an
lakemedel_10_funktion_Perorala_a
lakemedel_11_funktion_Perorala_a
lakemedel_12_funktion_Perorala_a
lakemedel_13_funktion_Perorala_a
lakemedel_14_funktion_Perorala_a
lakemedel_15_funktion_Perorala_a
lakemedel_16_funktion_Perorala_a
lakemedel_17_funktion_Perorala_a
lakemedel_18_funktion_Perorala_a
lakemedel_19_funktion_Perorala_a
lakemedel_20_funktion_Perorala_a;


do i = 1 to 20;
  if antihyp_d {i} = 0 or diuretic_d {i}=6 then antihyp_h=1;
  if insulin_d {i} = 2 or perorala_d {i}=1 then diab_treat=1;
end;
keep TWINNr diab_treat antihyp_h;

if antihyp_h=. then antihyp_h=0;
if diab_treat=. then diab_treat=0;

run;


data smoke1;
set temp2.Harmony_fritid;
if rokt_Ja_fortfarande=2 or rokt_Ja_men_ej_kontinuerligt=0 then smoke1=1; else smoke1=0;
keep TWINNR smoke1 ;
run;

* Smoke for subjects with missing smoke;
data smoke2;
set temp2.Harmony_fritid_informant;
if rokt_Ja_fortfarande=2 or rokt_Ja_men_ej_kontinuerligt=0 then smoke2=1; else smoke2=0;
keep TWINNR smoke2 ;
run;

proc sort data=smoke1; by twinnr; run;
proc sort data=smoke2; by twinnr; run;

data smoke;
merge smoke1 smoke2;
by TWINNR;
smoke_h=smoke1;
if smoke1 = . then smoke_h=smoke2;
keep TWINNR smoke_h;
run;



proc sort data=FHS; by TWINNR; run;
proc sort data=lab; by TWINNR; run;
proc sort data=diab; by TWINNR; run;
proc sort data=diab_treat_antihyp; by TWINNR; run;
proc sort data=smoke; by TWINNR; run;


** harmony;

data harmony;
merge FHS  lab diab diab_treat_antihyp smoke ;
by TWINNR;

*Indicator;
harmony=1;


* character to numeric;
tc_h=input(KOL,best12.);
hdl_h=input(HDL,best12.);
hba1=input(HBA1C,best12.);

* Calculate diabetes;
if diab_h = . and diab_treat = . then diab_h = .;
else if diab_h =1 or diab_treat=1 or HBA1C ge 6.5  then diab_h = 1;
else diab_h = 0;

* Keep only people with visit;
if check_h ne .;
drug_lipid_h=0;

* MISSING;
sbp_m=missing(sbp_h);
bmi_m=missing(bmi_h);
diab_m=missing(diab_h);
antihyp_m=missing(antihyp_h);
drug_lipid_m=missing(drug_lipid_h);
smoke_m=missing(smoke_h);
tc_m=missing(tc_h);
hdl_m=missing(hdl_h);

miss_h=sbp_m+bmi_m+diab_m+antihyp_m+smoke_m+tc_m+hdl_m+drug_lipid_m;

drop diab_treat KOL HDL HBA1C hba1 sbp_m bmi_m diab_m antihyp_m smoke_m tc_m hdl_m drug_lipid_m;
run;



******************************************************
************  MERGING ALL TOGETHER  ******************
******************************************************;

**IMPORT BESTZYGOSITY, SEX and PAIRID;
data infos;
set stradmin.BESTZYG;
keep TWINNR SEX BESTZYG PAIRID;
run;

**IMPORT BIRTHDATE;
data bdate;
set stradmin.PERSON_INFO;
birth_d=mdy(substr(BIRTHDATE ,5,2),substr(BIRTHDATE ,7,2), substr(BIRTHDATE ,1,4));
keep twinnr birth_d ;
format birth_d date9.;
run;


**IMPORT fam file so i keep only these subjects;
PROC IMPORT OUT= WORK.Fam 
            DATAFILE= "&path.\Data\GOSH\twins.+.qced.imputed.ceu.p2.fam" 
            DBMS=DLM REPLACE;
     DELIMITER=' '; 
     GETNAMES=NO;
RUN;

data fam;
set fam;
ind=0;
ind=substr(VAR2,1,6);
if ind not in (900,902);
TWINNR=VAR2;
keep TWINNR;
run;

*PRE - MERGING OF SATSA 1-2-3;


proc sort data=satsa_1; by TWINNR; run;
proc sort data=satsa_2; by TWINNR; run;
proc sort data=satsa_3; by TWINNR; run;


data satsa;
set satsa_3;
by twinnr;

* CREATE STUDY INDICATOR;
satsa=1;

study_s="satsa_3";

check_s=check_s_3;
bmi_s=bmi_s_3;
sbp_s=sbp_s_3;
smoke_s=smoke_s_3;
antihyp_s=antihyp_s_3;
drug_lipid_s=drug_lipid_s_3;
diab_s=diab_s_3;
tc_s=tc_s_3;
hdl_s=hdl_s_3;
miss_s=miss_s_3;

keep TWINNR check_s bmi_s sbp_s smoke_s antihyp_s tc_s hdl_s diab_s study_s satsa miss_s drug_lipid_s; 

format check_s date9.;
run;




*MERGE EVERYTHING;

proc sort data=fam; by TWINNR; run;
proc sort data=infos; by TWINNR; run;
proc sort data=bdate; by TWINNR; run;
proc sort data=gender; by TWINNR; run;
proc sort data=octo; by TWINNR; run;
proc sort data=satsa; by TWINNR; run;
proc sort data=harmony; by TWINNR; run;


data gosh.GOSH;
merge gender octo satsa harmony infos bdate fam(in=x);
by twinnr;
if x=1;

* CREATE STUDY INDICATOR;

if octo=. then octo=0;
if harmony=. then harmony=0;
if satsa=. then satsa=0;
if gender=. then gender=0;
count=octo+harmony+satsa+gender;

study="             ";

* FIX SUBJECTS IN MORE THAN ONE STUDY;

if count>1 then do;

  if miss_s = 0 then do;
      study="satsa";
      check_d=check_s;
  end;
  else if miss_o = 0 then do;
      study="octo";
      check_d=check_o;
  end;
  else if miss_g = 0 then do;
      study="gender";
      check_d=check_g;
  end;
  else if miss_o ne 0 and miss_s ne 0 and miss_g ne 0 then do;
      study="harmony";
      check_d=check_h;
  end;

  if miss_s = 0 then bmi=bmi_s;
  else if miss_o = 0 then bmi=bmi_o;
  else bmi=bmi_h;

  if miss_s = 0 then sbp=sbp_s;
  else if miss_o = 0 then sbp=sbp_o;
  else sbp=sbp_h;

  if miss_s = 0 then smoke=smoke_s;
  else if miss_o = 0 then smoke=smoke_o;
  else smoke=smoke_h;

  if miss_s = 0 then antihyp=antihyp_s;
  else if miss_o = 0 then antihyp=antihyp_o;
  else antihyp=antihyp_h;

  if miss_s = 0 then diab=diab_s;
  else if miss_o = 0 then diab=diab_o;
  else diab=diab_h;

  if miss_s = 0 then drug_lipid=drug_lipid_s;
  else if miss_o = 0 then drug_lipid=drug_lipid_o;
  else drug_lipid=drug_lipid_h;

  if miss_s = 0 then tc=tc_s;
  else if miss_o = 0 then tc=tc_o;
  else tc=tc_h;

  if miss_s = 0 then hdl=hdl_s;
  else if miss_o = 0 then hdl=hdl_o;
  else hdl=hdl_h;

 end;

* SUBJECTS IN ONLY 1 study;

if count=1 then do;

  if missing(check_g)=0  then check_d=check_g;
  if missing(check_o)=0  then check_d=check_o;
  if missing(check_s)=0  then check_d=check_s;
  if missing(check_h)=0  then check_d=check_h;

  if missing(bmi_g)=0  then bmi=bmi_g;
  if missing(bmi_o)=0  then bmi=bmi_o;
  if missing(bmi_s)=0  then bmi=bmi_s;
  if missing(bmi_h)=0  then bmi=bmi_h;

  if missing(sbp_g)=0  then sbp=sbp_g;
  if missing(sbp_o)=0  then sbp=sbp_o;
  if missing(sbp_s)=0  then sbp=sbp_s;
  if missing(sbp_h)=0  then sbp=sbp_h;

  if missing(smoke_g)=0  then smoke=smoke_g;
  if missing(smoke_o)=0  then smoke=smoke_o;
  if missing(smoke_s)=0  then smoke=smoke_s;
  if missing(smoke_h)=0  then smoke=smoke_h;

  if missing(antihyp_g)=0  then antihyp=antihyp_g;
  if missing(antihyp_o)=0  then antihyp=antihyp_o;
  if missing(antihyp_s)=0  then antihyp=antihyp_s;
  if missing(antihyp_h)=0  then antihyp=antihyp_h;

  if missing(drug_lipid_g)=0  then drug_lipid=drug_lipid_g;
  if missing(drug_lipid_o)=0  then drug_lipid=drug_lipid_o;
  if missing(drug_lipid_s)=0  then drug_lipid=drug_lipid_s;
  if missing(drug_lipid_h)=0  then drug_lipid=drug_lipid_h;

  if missing(tc_g)=0  then tc=tc_g;
  if missing(tc_o)=0  then tc=tc_o;
  if missing(tc_s)=0  then tc=tc_s;
  if missing(tc_h)=0  then tc=tc_h;

  if missing(apoa_g)=0  then hdl=apoa_g;
  if missing(hdl_o)=0  then hdl=hdl_o;
  if missing(hdl_s)=0  then hdl=hdl_s;
  if missing(hdl_h)=0  then hdl=hdl_h;

  if missing(diab_g)=0  then diab=diab_g;
  if missing(diab_o)=0  then diab=diab_o;
  if missing(diab_s)=0  then diab=diab_s;
  if missing(diab_h)=0  then diab=diab_h;

  if gender = 1 then study="gender";
  if octo = 1 then study="octo";
  if satsa = 1 then study="satsa";
  if harmony = 1 then study="harmony";
end;


keep TWINNR PAIRID SEX BESTZYG check_d birth_d bmi sbp smoke antihyp tc hdl diab study study_s drug_lipid; 

format check_d date9.;
run;


PROC EXPORT DATA= gosh.gosh
            OUTFILE= "&path.\Data\GOSH\Final\GOSH.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=YES;
RUN;



**********************
******** TESTS *******
**********************;

** Distribution analysis;
proc sort data=gosh.GOSH; by study; run;

proc freq data=gosh.GOSH;
by study;
table smoke antihyp diab ;
run;

proc means data=gosh.GOSH;
by study;
var bmi sbp tc hdl;
run;

** Number of subjects;
proc means data=gosh.GOSH n;
by study;
var TWINNR;
run;


** Associaiton;
proc phreg data=gosh.GOSH;
model (age_entry_chd,age_exit_chd)*incchd(0) = sex tc hdl antihyp sbp smoke bmi diab;
where chdd > check_d;
run;

** Compare with Erik's data;
PROC IMPORT OUT= WORK.erik 
            DATAFILE= "&path.\Data\GOSH\GOSH_metabochip_phenotype_master_dataset.dta" 
            DBMS=STATA REPLACE;

RUN;

* No duplicated twins;
proc freq data=erik noprint; table TWINNR / out=temp; run;

** Distribution analysis;
data erik; set erik; if cohort in (1,2,3) then cohort =1;
else cohort=cohort;
run;

proc sort data=erik; by cohort; run;

proc freq data=erik;
by cohort;
table smoking  bpmed selfrep_diab diabmed insul;
run;

