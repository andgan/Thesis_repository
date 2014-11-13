/*----------------------------------------------
Filename: merge_ULSAM_TWINGENE_GOSH_score_polygene_review.sas
Study: genscore
Author: Andrea Ganna
Date: 23FEB2013
Updated: 21JAN2013 - added export of only CHD SNPs for other project
Purpose: 1. Merge ULSAM TWINGENE and GOSH PHENOTYPES and check that TWINGENE and GOSH don't have same twins
         2. IMPORT THE GENOTYPES SCORES AND MERGE WITH THE FINAL DATASET.
         3. IMPORT POLYEGEN SCORE AND MERGE WITH THE FINAL DATASET.
Note: This is after review, so one new score has been added
-----------------------------------------------
Data used: twge.Twingene ulsam.Ulsam_70 gosh.Gosh_cvd gen.txt export_twingene.txt export_gosh.txt export_ulsam.txt
Data created: ULSAM_GOSH_TWGE_final.csv (for R) ULSAM_GOSH_TWGE_final.dta (for stata)
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/




**** ASSIGN LIBRARIES;

%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;
*%let path=E:\Phd_KI\Genscore;


libname twge "&path.\Data\Twingene\Final";
libname ulsam "&path.\Data\ULSAM\Final";
libname gosh "&path.\Data\GOSH\Final";


proc sort data=twge.Twingene; by TWINNR; run;
proc sort data=ulsam.Ulsam_70; by GWAS_ID; run;
proc sort data=gosh.Gosh_cvd; by GWAS_ID; run;


** FIRST CHECK IF TWINGENE AND GOSH HAVE TWINS IN COMMON (in GOSH TWINNR=GWAS_ID);
data gosh;
set  gosh.Gosh_cvd;
TWINNR=input(GWAS_ID, 12.);
run;

proc sort data=gosh; by twinnr; run;

* Doing like this it keeps variables values for gosh when a twin is both in twingene and in gosh; 
data twins;
merge  twge.Twingene (in=x)  gosh (in=y);
by 	TWINNR;
run;


** CREATE THE FINAL MERGED FILE;
data final;
set twins (drop=TWINNR) ulsam.Ulsam_70;
run;

*** THIS FOR EXPORTING FOR EXTRA ANALYSIS;

PROC IMPORT OUT= WORK.chd_snps 
            DATAFILE= "\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\From_plink\single_ana_chd.txt" 
            DBMS=DLM REPLACE;
     DELIMITER='00'x; 
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc sort data=chd_snps; by GWAS_ID; run;
proc sort data=final; by GWAS_ID; run;
data chd_final;
merge chd_snps final;
by 	GWAS_ID;
if check_d=. or birth_d=. then delete;
drop famhist;
run;


PROC EXPORT DATA= WORK.chd_final 
            OUTFILE= "&path.\Results\CHD_snps_ULSAM_GOSH_TWINGENE.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


** Check there are no duplicares;

proc freq data=final  noprint; table GWAS_ID / out=temp; run;


*** NOW IMPORT SCORES;

      data WORK.GEN                                     ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\From_plink\m_gen_review.txt'
    delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
         format GWAS_ID $20. ;
         format ALLCAT best12. ;
		 format ALLCAT_q best12. ;
         format ALLCAT_weights_cl_wtccc best12. ;
		 format ALLCAT_weights_cl_wtccc_q best12. ;
         format ALLCAT_bmi best12. ;
         format ALLCAT_bmi_q best12. ;
         format ALLCAT_bmi_cl_wtccc best12. ;
         format ALLCAT_bmi_cl_wtccc_q best12. ;
         format ALLCAT_chd best12. ;
         format ALLCAT_chd_q best12. ;
         format ALLCAT_chd_cl_wtccc best12. ;
         format ALLCAT_chd_cl_wtccc_q best12. ;
         format ALLCAT_hdl best12. ;
         format ALLCAT_hdl_q best12. ;
         format ALLCAT_hdl_cl_wtccc best12. ;
         format ALLCAT_hdl_cl_wtccc_q best12. ;
         format ALLCAT_sbp best12. ;
         format ALLCAT_sbp_q best12. ;
         format ALLCAT_sbp_cl_wtccc best12. ;
         format ALLCAT_sbp_cl_wtccc_q best12. ;
         format ALLCAT_smoke best12. ;
         format ALLCAT_smoke_q best12. ;
         format ALLCAT_smoke_cl_wtccc best12. ;
         format ALLCAT_smoke_cl_wtccc_q best12. ;
         format ALLCAT_t2d best12. ;
         format ALLCAT_t2d_q best12. ;
         format ALLCAT_t2d_cl_wtccc best12. ;
         format ALLCAT_t2d_cl_wtccc_q best12. ;
         format ALLCAT_tc best12. ;
         format ALLCAT_tc_q best12. ;
         format ALLCAT_tc_cl_wtccc best12. ;
         format ALLCAT_tc_cl_wtccc_q best12. ;
		 format ALLCAT_fhs best12. ;
         format ALLCAT_fhs_q best12. ;
		 format ALLCAT_fhs_cl_wtccc best12. ;
         format ALLCAT_fhs_cl_wtccc_q best12. ;
         format ALLCAT_cl_wtccc best12. ;
         format ALLCAT_cl_wtccc_q best12. ;
         format ALLCAT_chd_weights_cl_wtccc best12. ;
		 format ALLCAT_chd_weights_cl_wtccc_q best12. ;
         format STAR_chd best12. ;
         format STAR_chd_q best12. ;
		 format allcat_no_chd_weights_cl_wtccc best12. ;
         format allcat_no_chd_weights_cl_wtccc_q best12. ;
		 
        
                  
      input
                  GWAS_ID $
                  ALLCAT
				  ALLCAT_q
                  ALLCAT_weights_cl_wtccc
				  ALLCAT_weights_cl_wtccc_q
                  ALLCAT_bmi
                  ALLCAT_bmi_q
                  ALLCAT_bmi_cl_wtccc
                  ALLCAT_bmi_cl_wtccc_q
                  ALLCAT_chd
                  ALLCAT_chd_q
                  ALLCAT_chd_cl_wtccc
                  ALLCAT_chd_cl_wtccc_q
                  ALLCAT_hdl
                  ALLCAT_hdl_q
                  ALLCAT_hdl_cl_wtccc
                  ALLCAT_hdl_cl_wtccc_q
                  ALLCAT_sbp
                  ALLCAT_sbp_q
                  ALLCAT_sbp_cl_wtccc
                  ALLCAT_sbp_cl_wtccc_q
                  ALLCAT_smoke
                  ALLCAT_smoke_q
                  ALLCAT_smoke_cl_wtccc
                  ALLCAT_smoke_cl_wtccc_q
                  ALLCAT_t2d
                  ALLCAT_t2d_q
                  ALLCAT_t2d_cl_wtccc
                  ALLCAT_t2d_cl_wtccc_q
                  ALLCAT_tc
                  ALLCAT_tc_q
                  ALLCAT_tc_cl_wtccc
                  ALLCAT_tc_cl_wtccc_q
                  ALLCAT_fhs
                  ALLCAT_fhs_q
                  ALLCAT_fhs_cl_wtccc
                  ALLCAT_fhs_cl_wtccc_q
                  ALLCAT_cl_wtccc
                  ALLCAT_cl_wtccc_q
                  ALLCAT_chd_weights_cl_wtccc 
		          ALLCAT_chd_weights_cl_wtccc_q 
                  STAR_chd
                  STAR_chd_q
				  allcat_no_chd_weights_cl_wtccc
				  allcat_no_chd_weights_cl_wtccc_q
      ;
      run;

*** IMPORT POLYGENE SCORES;

PROC IMPORT OUT= WORK.p_tw 
            DATAFILE= "&path.\Results\polygene\export_twingene.txt" 
            DBMS=DLM REPLACE;
     DELIMITER=' '; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data p_tw; set p_tw; GWAS_ID=var18; drop var18; run;

PROC IMPORT OUT= WORK.p_go 
            DATAFILE= "&path.\Results\polygene\export_gosh.txt" 
            DBMS=DLM REPLACE;
     DELIMITER=' '; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data p_go; set p_go; GWAS_ID=compress(put(var18, 20.)); drop var18; run;

PROC IMPORT OUT= WORK.p_ul 
            DATAFILE= "&path.\Results\polygene\export_ulsam.txt" 
            DBMS=DLM REPLACE;
     DELIMITER=' '; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data p_ul; set p_ul; GWAS_ID="UL"||compress(put(var18, 18.));drop var18; run;

data p;
set p_go p_tw p_ul;
run;

** NOW MERGE;

proc sort data=gen; by GWAS_ID; run;
proc sort data=final; by GWAS_ID; run;
proc sort data=p; by GWAS_ID; run;

data final_gen;
merge final gen (in=x) p;
by GWAS_ID;
if x=1;
* This subjects were the duplicated subjects in both twingene and gosh and are now deleted from ge;
if GWAS_FID ne "";

* Subjects without birth date are deleted;
if birth_d = . then delete;

* Determine a CVD event (first between CHD and IS);
if chdd <= isd then do;
	cvdd=chdd;
	cvd=chd;
end;
else if chdd > isd then do;
	cvdd=isd;
	cvd=is;
end;

* Incident CVD;
if cvdd > check_d and cvd=1 then inccvd=1; else inccvd=0;

* Age entry and exit CVD;

age_entry_cvd = age;
age_exit_cvd = 	(cvdd - birth_d) / 365.25;

* Outliers;
if age>110 or age<20 then age=.;
if sbp>250 or sbp<70 then sbp=.;
if bmi>49.2 or bmi<5.8 then bmi=.;

if study="" then study="gender";

format cvdd date9.;
run;


** Export for STATA analysis;

PROC EXPORT DATA= WORK.FINAL_GEN 
            OUTFILE= "&path.\Results\ULSAM_GOSH_TWGE_final.dta" 
            DBMS=DTA REPLACE;
RUN;

** Export for R;

PROC EXPORT DATA= WORK.FINAL_GEN 
            OUTFILE= "&path.\Results\ULSAM_GOSH_TWGE_final_review.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

** Export only phenotypes;

data only_use;
set FINAL_GEN;
keep GWAS_ID GWAS_FID age_entry_chd age_exit_chd SEX tc hdl smoke antihyp sbp diab bmi study incchd check_d chdd;
run; 

PROC EXPORT DATA= WORK.only_use 
            OUTFILE= "&path.\Results\ULSAM_GOSH_TWGE_final_only_use.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


** DESCRIPTIVE STATISTICS;

proc freq data=final_gen;
table smoke antihyp diab ;
run;

proc means data=final_gen;
var bmi sbp tc hdl;
run;

** Association;

proc phreg data=final_gen;
model (age_entry_chd,age_exit_chd)*incchd(0) = sex tc hdl antihyp sbp smoke bmi diab famhist;
where chdd > check_d;
run;

