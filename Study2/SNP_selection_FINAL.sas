/*----------------------------------------------
Filename: SNP_selection_FINAL.sas
Study: Genscore
Author: Andrea Ganna
Date: 14JAN2013
Updated: 23FEN2013, becuase of a reviewer I export also the allcat score without CHD-associated SNPs
Purpose: Import all SNPs from catalog; clean and format data
Note:
-----------------------------------------------
Data used: NHGRI_140113.txt, NHGRI_add_DIAB.csv, NHGRI_add_CAD.csv
Data created: LOTS OF SCORES (Check index)
-----------------------------------------------
OP: SAS 9.2
-----------------------------------------------*/

***********************	 
******** INDEX ********
***********************


1. IMPORT DATA FROM THE CATALOG 
2. CHECK RACE AND SAMPLE SIZE
3. Add CHD for sensitivity
4. MERGE, CREATE GROUPS AND SOLVE MINOR PROBLEMS 
5. CHECK: MISSING RISK ALLELE OR ALLELE FREQUENCY
6. ADD PROXIES FOR SNPs not FOUND in HAPMAP relase 27 and in TWINGENE GWAS
7. ANNOTATION 
-------- PGM BREAK -------
8. CHECK: ANNOTATION PROBLEMS
9. FLIP - DEFINITIVE DATASET
10. DUPLICATED SNPs within  TRAITS
11. SNPs quality metrics
12. EXPORT SNPLIST

* Specific scores
e1. EXPORT: CATALOG SCORE --> ALLCAT
e2.	EXPORT: STATE-OF-ART SCORE --> STAR
e3. EXPORT: FHS specific trait from catalog
e4. EXPORT: FHS specific trait from state-of-art papers
;


%let path=\\.psf\Home\Documents\Work\Phd_KI\Genscore;

* Settings libnames;

libname sel "&path.\Data\SNPs_selection";

*****************************************************
1.c. IMPORT DATA FROM THE CATALOG (UPDATED 10OCT2011)* 
*****************************************************;
     data WORK.NHGRI                                    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\SNPs_selection\NHGRI_100712.txt' delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat Date_Added_to_Catalog mmddyy10. ;
         informat PUBMEDID best32. ;
         informat First_Author $15. ;
         informat Date mmddyy10. ;
         informat Journal $23. ;
         informat Link $43. ;
         informat Study $243. ;
         informat Disease_Trait $80. ;
         informat Initial_Sample_Size $216. ;
         informat Replication_Sample_Size $226. ;
         informat Region $9. ;
         informat Chr_id best32. ;
         informat Chr_pos best32. ;
         informat Reported_Gene_s_ $55. ;
         informat Mapped_gene $55. ;
         informat Upstream_gene_id best32. ;
         informat Downstream_gene_id best32. ;
         informat Snp_gene_ids $17. ;
         informat Upstream_gene_distance best32. ;
         informat Downstream_gene_distance best32. ;
         informat Strongest_SNP_Risk_Allele $15. ;
         informat SNPs $13. ;
         informat Merged best32. ;
         informat Snp_id_current $13. ;
         informat Context $21. ;
         informat Intergenic best32. ;
         informat Risk_Allele_Frequency $6. ;
         informat p_Value $6. ;
         informat Pvalue_mlog $6. ;
         informat p_Value__text_ $27. ;
         informat OR_or_beta $5. ;
         informat _5__CI__text_ $40. ;
         informat Platform__SNPs_passing_QC_ $55. ;
         informat CNV $1. ;
         format Date_Added_to_Catalog mmddyy10. ;
         format PUBMEDID best12. ;
         format First_Author $15. ;
         format Date mmddyy10. ;
         format Journal $23. ;
         format Link $43. ;
         format Study $243. ;
         format Disease_Trait $80. ;
         format Initial_Sample_Size $215. ;
         format Replication_Sample_Size $226. ;
         format Region $9. ;
         format Chr_id best32. ;
         format Chr_pos best32. ;
         format Reported_Gene_s_ $55. ;
         format Mapped_gene $55. ;
         format Upstream_gene_id best32. ;
         format Downstream_gene_id best32. ;
         format Snp_gene_ids $17. ;
         format Upstream_gene_distance best32. ;
         format Downstream_gene_distance best32. ;
         format Strongest_SNP_Risk_Allele $15. ;
         format SNPs $13. ;
         format Merged best32. ;
         format Snp_id_current $13. ;
         format Context $21. ;
         format Intergenic best32. ;
         format Risk_Allele_Frequency $6. ;
         format p_Value $6. ;
         format Pvalue_mlog $6. ;
         format p_Value__text_ $27. ;
         format OR_or_beta $5. ;
         format _5__CI__text_ $40. ;
         format Platform__SNPs_passing_QC_ $55. ;
         format CNV $1. ;
      input
                  Date_Added_to_Catalog
                  PUBMEDID
                  First_Author $
                  Date
                  Journal $
                  Link $
                  Study $
                  Disease_Trait $
                  Initial_Sample_Size $
                  Replication_Sample_Size $
                  Region $
                  Chr_id
                  Chr_pos
                  Reported_Gene_s_ $
                  Mapped_gene $
                  Upstream_gene_id
                  Downstream_gene_id
                  Snp_gene_ids $
                  Upstream_gene_distance
                  Downstream_gene_distance
                  Strongest_SNP_Risk_Allele $
                  SNPs $
                  Merged
                  Snp_id_current $
                  Context $
                  Intergenic
                  Risk_Allele_Frequency $
                  p_Value $
                  Pvalue_mlog $
                  p_Value__text_ $
                  OR_or_beta $
                 _5__CI__text_ $
                  Platform__SNPs_passing_QC_ $
                  CNV $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;


data WORK.NHGRI2                                    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\SNPs_selection\NHGRI_140113.txt' delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat Date_Added_to_Catalog mmddyy10. ;
         informat PUBMEDID best32. ;
         informat First_Author $15. ;
         informat Date mmddyy10. ;
         informat Journal $23. ;
         informat Link $43. ;
         informat Study $243. ;
         informat Disease_Trait $80. ;
         informat Initial_Sample_Size $216. ;
         informat Replication_Sample_Size $226. ;
         informat Region $9. ;
         informat Chr_id best32. ;
         informat Chr_pos best32. ;
         informat Reported_Gene_s_ $55. ;
         informat Mapped_gene $55. ;
         informat Upstream_gene_id best32. ;
         informat Downstream_gene_id best32. ;
         informat Snp_gene_ids $17. ;
         informat Upstream_gene_distance best32. ;
         informat Downstream_gene_distance best32. ;
         informat Strongest_SNP_Risk_Allele $15. ;
         informat SNPs $13. ;
         informat Merged best32. ;
         informat Snp_id_current $13. ;
         informat Context $21. ;
         informat Intergenic best32. ;
         informat Risk_Allele_Frequency $6. ;
         informat p_Value $6. ;
         informat Pvalue_mlog $6. ;
         informat p_Value__text_ $27. ;
         informat OR_or_beta $5. ;
         informat _5__CI__text_ $40. ;
         informat Platform__SNPs_passing_QC_ $55. ;
         informat CNV $1. ;
         format Date_Added_to_Catalog mmddyy10. ;
         format PUBMEDID best12. ;
         format First_Author $15. ;
         format Date mmddyy10. ;
         format Journal $23. ;
         format Link $43. ;
         format Study $243. ;
         format Disease_Trait $80. ;
         format Initial_Sample_Size $215. ;
         format Replication_Sample_Size $226. ;
         format Region $9. ;
         format Chr_id best32. ;
         format Chr_pos best32. ;
         format Reported_Gene_s_ $55. ;
         format Mapped_gene $55. ;
         format Upstream_gene_id best32. ;
         format Downstream_gene_id best32. ;
         format Snp_gene_ids $17. ;
         format Upstream_gene_distance best32. ;
         format Downstream_gene_distance best32. ;
         format Strongest_SNP_Risk_Allele $15. ;
         format SNPs $13. ;
         format Merged best32. ;
         format Snp_id_current $13. ;
         format Context $21. ;
         format Intergenic best32. ;
         format Risk_Allele_Frequency $6. ;
         format p_Value $6. ;
         format Pvalue_mlog $6. ;
         format p_Value__text_ $27. ;
         format OR_or_beta $5. ;
         format _5__CI__text_ $40. ;
         format Platform__SNPs_passing_QC_ $55. ;
         format CNV $1. ;
      input
                  Date_Added_to_Catalog
                  PUBMEDID
                  First_Author $
                  Date
                  Journal $
                  Link $
                  Study $
                  Disease_Trait $
                  Initial_Sample_Size $
                  Replication_Sample_Size $
                  Region $
                  Chr_id
                  Chr_pos
                  Reported_Gene_s_ $
                  Mapped_gene $
                  Upstream_gene_id
                  Downstream_gene_id
                  Snp_gene_ids $
                  Upstream_gene_distance
                  Downstream_gene_distance
                  Strongest_SNP_Risk_Allele $
                  SNPs $
                  Merged
                  Snp_id_current $
                  Context $
                  Intergenic
                  Risk_Allele_Frequency $
                  p_Value $
                  Pvalue_mlog $
                  p_Value__text_ $
                  OR_or_beta $
                 _5__CI__text_ $
                  Platform__SNPs_passing_QC_ $
                  CNV $
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;

* Check differences with old catalog;
proc freq data=NHGRI noprint; table Disease_Trait / out=test; run;
proc freq data=NHGRI2 noprint; table Disease_Trait / out=test2; run;

data tests;
merge test (in=x) test2 (in=y);
by Disease_Trait;
if x=0 and y=1;
run;


  data WORK.NHGRI_P                                 ;
 infile '\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\SNPs_selection\NHGRI_add_CAD.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat First_Author $8. ;
         informat Journal $23. ;
         informat Study $243. ;
         informat Disease_Trait $80. ;
         informat Initial_Sample_Size $215. ;
         informat Replication_Sample_Size $226. ;
         informat Reported_Gene_s_ $55. ;
         informat Strongest_SNP_Risk_Allele $15. ;
         informat SNPs $13. ;
         informat Risk_Allele_Frequency $6. ;
         informat p_Value $6. ;
         informat OR_or_beta $5. ;
         informat _5__CI__text_ $14. ;
         informat PUBMEDID best32. ;
         format First_Author $8. ;
         format Journal $23. ;
         format Study $243. ;
         format Disease_Trait $80. ;
         format Initial_Sample_Size $215. ;
         format Replication_Sample_Size $226. ;
         format Reported_Gene_s_ $55. ;
         format Strongest_SNP_Risk_Allele $15. ;
         format SNPs $13. ;
         format Risk_Allele_Frequency $6. ;
         format p_Value $6. ;
         format OR_or_beta $5. ;
         format _5__CI__text_ $14. ;
         format PUBMEDID best12. ;
      input
                  First_Author $
                  Journal $
                  Study $
                  Disease_Trait $
                  Initial_Sample_Size $
                  Replication_Sample_Size $
                  Reported_Gene_s_ $
                  Strongest_SNP_Risk_Allele $
                  SNPs $
                  Risk_Allele_Frequency	$
                  p_Value $
                  OR_or_beta $
                  _5__CI__text_ $
                  PUBMEDID
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;



  data WORK.NHGRI_P2                                 ;
 infile '\\psf\Home\Documents\Work\Phd_KI\Genscore\Data\SNPs_selection\NHGRI_add_DIAB.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat First_Author $8. ;
         informat Journal $23. ;
         informat Study $243. ;
         informat Disease_Trait $80. ;
         informat Initial_Sample_Size $215. ;
         informat Replication_Sample_Size $226. ;
         informat Reported_Gene_s_ $55. ;
         informat Strongest_SNP_Risk_Allele $15. ;
         informat SNPs $13. ;
         informat Risk_Allele_Frequency $6. ;
         informat p_Value $6. ;
         informat OR_or_beta $5. ;
         informat _5__CI__text_ $14. ;
         informat PUBMEDID best32. ;
         format First_Author $8. ;
         format Journal $23. ;
         format Study $243. ;
         format Disease_Trait $80. ;
         format Initial_Sample_Size $215. ;
         format Replication_Sample_Size $226. ;
         format Reported_Gene_s_ $55. ;
         format Strongest_SNP_Risk_Allele $15. ;
         format SNPs $13. ;
         format Risk_Allele_Frequency $6. ;
         format p_Value $6. ;
         format OR_or_beta $5. ;
         format _5__CI__text_ $14. ;
         format PUBMEDID best12. ;
      input
                  First_Author $
                  Journal $
                  Study $
                  Disease_Trait $
                  Initial_Sample_Size $
                  Replication_Sample_Size $
                  Reported_Gene_s_ $
                  Strongest_SNP_Risk_Allele $
                  SNPs $
                  Risk_Allele_Frequency	$
                  p_Value $
                  OR_or_beta $
                  _5__CI__text_ $
                  PUBMEDID
      ;
      if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
      run;



data NHGRI_;
set NHGRI2 NHGRI_p NHGRI_p2;
run;




** Check traits to include in the analysis;
proc freq data=NHGRI_ noprint;
table Disease_Trait/ out=NHGRI_trait;
run;


proc format;

value traitf 1= "abdominal aortic aneurysm"
2="adiponectin levels"
3="adiposity"
4="thoracic aortic aneurysms and dissections"
5="aortic root size"
6="arterial stiffness"
7="atrial fibrillation"
8="atrial fibrillation/atrial flutter"
9="atrioventricular conduction"
11="blood pressure"
12="body mass (lean)"
13="body mass index"
14="body mass index and fat mass"
15="c-reactive protein"
16="cardiac structure and function"
17="cardiovascular disease risk factors"
18="cholesterol"
19="cholesterol, total"
20="coronary artery calcification"
21="coronary heart disease"
22="coronary spasm"
23="diabetes (incident)"
24="diabetes related insulin traits"
25="diastolic blood pressure"
26="echocardiographic traits"
27="electrocardiographic conduction measures"
28="electrocardiographic traits"
29="endothelial function traits"
30="fasting glucose-related traits"
31="fasting insulin-related traits"
32="fasting plasma glucose"
33="fibrinogen"
34="glycated hemoglobin levels"
35="hdl cholesterol"
36="heart failure"
37="heart rate variability traits"
38="hemostatic factors and hematological phenotypes"
39="hypertension"
40="hypertension (young onset)"
41="hypertriglyceridemia"
42="insulin resistance/response"
44="interleukin-18 levels"
45="ischemic stroke"
46="ldl cholesterol"
47="left ventricular mass"
48="proinsulin levels"
49="lipoprotein-associated phospholipase A2 activity and mass"
50="major cvd"
51="metabolic syndrome"
52="myocardial infarction"
53="myocardial infarction (early onset)"
54="obesity"
55="obesity (early onset extreme)"
56="obesity (extreme)"
57="obesity-related traits"
58="pr interval"
59="peripheral artery disease"
60="qt interval"
61="rr interval (heart rate)"
62="resting heart rate"
63="smoking behavior"
64="stroke"
65="subarachnoid aneurysmal hemorrhage"
66="subclinical atherosclerosis traits"
67="subclinical brain infarct"
68="sudden cardiac arrest"
69="systolic blood pressure"
70="triglycerides"
71="two-hour glucose challenge"
72="type 1 diabetes"
73="type 2 diabetes"
74="type 2 diabetes and 6 quantitative traits"
75="type 2 diabetes and other traits"
76="ventricular conduction"
77="ventricular fibrillation"
78="waist circumference"
79="waist circumference and related phenotypes"
80="waist-hip ratio"
81="weight"
83="cardiac hypertrophy"
84="dilated cardiomyopathy"
85="hdl cholesterol - triglycerides (hdlg-tg)"
86="metabolic syndrome (bivariate traits)"
87="peripartum cardiomyopathy"
88="sick sinus syndrome"
89="triglycerides-blood pressure (tg-bp)"
90="waist circumference - triglycerides (wc-tg)"
91="natriuretic peptide levels"
92="coronary restenosis"
93="carotid intima media thickness"
94="metabolic traits"
95="ankle-brachial index"
96="aortic stiffness"
97="body mass index and cholesterol (psychopharmacological treatment)"
98="cardiac repolarization"
99="creatinine levels"
100="diabetes (gestational)"
101="hypertension risk in short sleep duration"
102="inflammatory biomarkers"
103="insulin-related traits"
104="lipid levels in hepatitis c treatment"
105="lipid metabolism phenotypes"
106="lipid traits"
107="lp (a) levels"
108="matrix metalloproteinase levels"
109="obesity and blood pressure"
110="polyunsaturated fatty acid levels"
111="postoperative ventricular dysfunction"
112="preeclampsia"
113="type 2 diabetes and gout"
114="urate levels"
115="uric acid levels"
116="vwf and fviii levels"
117="stroke (ischemic)"
118="pericardial fat"
119="plasminogen activator inhibitor type 1 levels (pai-1)"
120="subcutaneous adipose tissue"
121="visceral adipose tissue adjusted for bmi"
122="visceral adipose tissue/subcutaneous adipose tissue ratio"
123="visceral fat";


run;

data step1;
set NHGRI_;
Disease_Trait = lowcase(Disease_Trait);
if compare(Disease_Trait , "abdominal aortic aneurysm")=0 then trait=1;
if compare(Disease_Trait , "adiponectin levels")=0  then trait=2;
if compare(Disease_Trait , "adiposity")=0  then trait=3;
if compare(Disease_Trait , "thoracic aortic aneurysms and dissections")=0  then trait=4;
if compare(Disease_Trait , "aortic root size")=0  then trait=5;
if compare(Disease_Trait , "arterial stiffness")=0  then trait=6;
if compare(Disease_Trait , "atrial fibrillation")=0  then trait=7;
if compare(Disease_Trait , "atrial fibrillation/atrial flutter")=0  then trait=8;
if compare(Disease_Trait , "atrioventricular conduction")=0  then trait=9;
if compare(Disease_Trait , "blood pressure")=0  then trait=11;
if compare(Disease_Trait , "body mass (lean)")=0  then trait=12;
if compare(Disease_Trait , "body mass index")=0  then trait=13;
if compare(Disease_Trait , "body mass index and fat mass")=0  then trait=14;
if compare(Disease_Trait , "c-reactive protein")=0  then trait=15;
if compare(Disease_Trait , "cardiac structure and function")=0  then trait=16;
if compare(Disease_Trait , "cardiovascular disease risk factors")=0  then trait=17;
if compare(Disease_Trait , "cholesterol")=0  then trait=18;
if compare(Disease_Trait , "cholesterol, total")=0  then trait=19;
if compare(Disease_Trait , "coronary artery calcification")=0  then trait=20;
if compare(Disease_Trait , "coronary heart disease")=0  then trait=21;
if compare(Disease_Trait , "coronary spasm")=0  then trait=22;
if compare(Disease_Trait , "diabetes (incident)")=0  then trait=23;
if compare(Disease_Trait , "diabetes related insulin traits")=0  then trait=24;
if compare(Disease_Trait , "diastolic blood pressure")=0  then trait=25;
if compare(Disease_Trait , "echocardiographic traits")=0  then trait=26;
if compare(Disease_Trait , "electrocardiographic conduction measures")=0  then trait=27;
if compare(Disease_Trait , "electrocardiographic traits")=0  then trait=28;
if compare(Disease_Trait , "endothelial function traits")=0  then trait=29;
if compare(Disease_Trait , "fasting glucose-related traits")=0  then trait=30;
if compare(Disease_Trait , "fasting insulin-related traits")=0  then trait=31;
if compare(Disease_Trait , "fasting plasma glucose")=0  then trait=32;
if compare(Disease_Trait , "fibrinogen")=0  then trait=33;
if compare(Disease_Trait , "glycated hemoglobin levels")=0  then trait=34;
if compare(Disease_Trait , "hdl cholesterol")=0  then trait=35;
if compare(Disease_Trait , "heart failure")=0  then trait=36;
if compare(Disease_Trait , "heart rate variability traits")=0  then trait=37;
if compare(Disease_Trait , "hemostatic factors and hematological phenotypes")=0  then trait=38;
if compare(Disease_Trait , "hypertension")=0  then trait=39;
if compare(Disease_Trait , "hypertension (young onset)")=0  then trait=40;
if compare(Disease_Trait , "hypertriglyceridemia")=0  then trait=41;
if compare(Disease_Trait , "insulin resistance/response")=0  then trait=42;
if compare(Disease_Trait , "interleukin-18 levels")=0  then trait=44;
if compare(Disease_Trait , "ischemic stroke")=0  then trait=45;
if compare(Disease_Trait , "ldl cholesterol")=0  then trait=46;
if compare(Disease_Trait , "left ventricular mass")=0  then trait=47;
if compare(Disease_Trait , "proinsulin levels")=0  then trait=48;
if compare(Disease_Trait , "lipoprotein-associated phospholipase a2 activity and mass")=0  then trait=49;
if compare(Disease_Trait , "major cvd")=0  then trait=50;
if compare(Disease_Trait , "metabolic syndrome")=0  then trait=51;
if compare(Disease_Trait , "myocardial infarction")=0  then trait=52;
if compare(Disease_Trait , "myocardial infarction (early onset)")=0  then trait=53;
if compare(Disease_Trait , "obesity")=0  then trait=54;
if compare(Disease_Trait , "obesity (early onset extreme)")=0  then trait=55;
if compare(Disease_Trait , "obesity (extreme)")=0  then trait=56;
if compare(Disease_Trait , "obesity-related traits")=0  then trait=57;
if compare(Disease_Trait , "pr interval")=0  then trait=58;
if compare(Disease_Trait , "peripheral artery disease")=0  then trait=59;
if compare(Disease_Trait , "qt interval")=0  then trait=60;
if compare(Disease_Trait , "rr interval (heart rate)")=0  then trait=61;
if compare(Disease_Trait , "resting heart rate")=0  then trait=62;
if compare(Disease_Trait , "smoking behavior")=0  then trait=63;
if compare(Disease_Trait , "stroke")=0  then trait=64;
if compare(Disease_Trait , "subarachnoid aneurysmal hemorrhage")=0  then trait=65;
if compare(Disease_Trait , "subclinical atherosclerosis traits (other)")=0  then trait=66;
if compare(Disease_Trait , "subclinical brain infarct")=0  then trait=67;
if compare(Disease_Trait , "sudden cardiac arrest")=0  then trait=68;
if compare(Disease_Trait , "systolic blood pressure")=0  then trait=69;
if compare(Disease_Trait , "triglycerides")=0  then trait=70;
if compare(Disease_Trait , "two-hour glucose challenge")=0  then trait=71;
if compare(Disease_Trait , "type 1 diabetes")=0  then trait=72;
if compare(Disease_Trait , "type 2 diabetes")=0  then trait=73;
if compare(Disease_Trait , "type 2 diabetes and 6 quantitative traits")=0  then trait=74;
if compare(Disease_Trait , "type 2 diabetes and other traits")=0  then trait=75;
if compare(Disease_Trait , "ventricular conduction")=0  then trait=76;
if compare(Disease_Trait , "ventricular fibrillation")=0  then trait=77;
if compare(Disease_Trait , "waist circumference")=0  then trait=78;
if compare(Disease_Trait , "waist circumference and related phenotypes")=0  then trait=79;
if compare(Disease_Trait , "waist-hip ratio")=0  then trait=80;
if compare(Disease_Trait , "weight")=0  then trait=81;
if compare(Disease_Trait , "cardiac hypertrophy")=0  then trait=83;
if compare(Disease_Trait , "dilated cardiomyopathy")=0  then trait=84;
if compare(Disease_Trait , "hdl cholesterol - triglycerides (hdlc-tg)")=0  then trait=85;
if compare(Disease_Trait , "metabolic syndrome (bivariate traits)")=0  then trait=86;
if compare(Disease_Trait , "peripartum cardiomyopathy")=0  then trait=87;
if compare(Disease_Trait , "sick sinus syndrome")=0  then trait=88;
if compare(Disease_Trait , "triglycerides-blood pressure (tg-bp)")=0  then trait=89;
if compare(Disease_Trait , "waist circumference - triglycerides (wc-tg)")=0  then trait=90;
if compare(Disease_Trait , "natriuretic peptide levels")=0  then trait=91;
if compare(Disease_Trait , "coronary restenosis")=0  then trait=92;
if compare(Disease_Trait , "carotid intima media thickness")=0  then trait=93;
if compare(Disease_Trait , "metabolic traits")=0  then trait=94;
if compare(Disease_Trait , "ankle-brachial index")=0  then trait=95;
if compare(Disease_Trait , "aortic stiffness")=0  then trait=96;
if compare(Disease_Trait , "body mass index and cholesterol (psychopharmacological treatment)")=0  then trait=97;
if compare(Disease_Trait , "cardiac repolarization")=0  then trait=98;
if compare(Disease_Trait , "creatinine levels")=0  then trait=99;
if compare(Disease_Trait , "diabetes (gestational)")=0  then trait=100;
if compare(Disease_Trait , "hypertension risk in short sleep duration")=0  then trait=101;
if compare(Disease_Trait , "inflammatory biomarkers")=0  then trait=102;
if compare(Disease_Trait , "insulin-related traits")=0  then trait=103;
if compare(Disease_Trait , "lipid levels in hepatitis c treatment")=0  then trait=104;
if compare(Disease_Trait , "lipid metabolism phenotypes")=0  then trait=105;
if compare(Disease_Trait , "lipid traits")=0  then trait=106;
if compare(Disease_Trait , "lp (a) levels")=0  then trait=107;
if compare(Disease_Trait , "matrix metalloproteinase levels")=0  then trait=108;
if compare(Disease_Trait , "obesity and blood pressure")=0  then trait=109;
if compare(Disease_Trait , "polyunsaturated fatty acid levels")=0  then trait=110;
if compare(Disease_Trait , "postoperative ventricular dysfunction")=0  then trait=111;
if compare(Disease_Trait , "preeclampsia")=0  then trait=112;
if compare(Disease_Trait , "type 2 diabetes and gout")=0  then trait=113;
if compare(Disease_Trait , "urate levels")=0  then trait=114;
if compare(Disease_Trait , "uric acid levels")=0  then trait=115;
if compare(Disease_Trait , "vwf and fviii levels")=0  then trait=116;
if compare(Disease_Trait , "stroke (ischemic)")=0  then trait=117;
if compare(Disease_Trait , "pericardial fat")=0  then trait=118;
if compare(Disease_Trait , "plasminogen activator inhibitor type 1 levels (pai-1)")=0  then trait=119;
if compare(Disease_Trait , "subcutaneous adipose tissue")=0  then trait=120;
if compare(Disease_Trait , "visceral adipose tissue adjusted for bmi")=0  then trait=121;
if compare(Disease_Trait , "visceral adipose tissue/subcutaneous adipose tissue ratio")=0  then trait=122;
if compare(Disease_Trait , "visceral fat")=0  then trait=123;


format trait traitf.;

if trait=. then delete;
run;


* Count number of SNPs;
proc freq data=step1 noprint;
table trait/ out=temp;
run;

proc freq data=step1 noprint;
table SNPs/ out=temp;
run;

proc freq data=step1 noprint;
table Study/ out=temp;
run;

*********************************
2. CHECK RACE AND SAMPLE SIZE	* 
*********************************;


proc freq data=step1 noprint;
table Initial_Sample_Size  / out=temp;
run;


*Exclude race and pvalues < 5E-8;

data step2;
set step1;

*Define that these are catalog SNPs;
NHGRI=1;


*Manual correction of sample size numbers;
if compare(Initial_Sample_Size,"100 > 445ms100 < 386ms")=0 then Initial_Sample_Size="100;100";
if compare(Initial_Sample_Size,"1,341-1,345  individuals, depending on measure(Framingham)")=0 then Initial_Sample_Size="1345";
if compare(Initial_Sample_Size,"548-1,175 individuals, depending on measure(Framingham)")=0 then Initial_Sample_Size="1175";
if compare(Initial_Sample_Size,"4,611 individuals (2,617 smokers)")=0 then Initial_Sample_Size="4611";
if compare(Initial_Sample_Size,"644-1,327 individuals, depending on measure(Framingham)")=0 then Initial_Sample_Size="1327";
if compare(Initial_Sample_Size,"673-984 individuals, depending on measure (Framingham)")=0 then Initial_Sample_Size="984";
if compare(Initial_Sample_Size,"(see Samani 2007)")=0 then Initial_Sample_Size="875 cases,1,644 controls";
if compare(Initial_Sample_Size,"(see Todd 2007)")=0 then Initial_Sample_Size="2997 trios, 4,000 cases, 5,000 controls";
if compare(Initial_Sample_Size,"(see Zeggini 2007)")=0 then Initial_Sample_Size="24,194 cases,55,598 controls";
if compare(Initial_Sample_Size,"2,033 individuals in 519 families; 1,461 twins (1/pair selected randomly)")=0 then Initial_Sample_Size="2,033 individuals in 519 families; 1,461 twins";
if compare(Initial_Sample_Size,"2,269 individuals in 644 families")=0 then Initial_Sample_Size="2,269 individuals";
if compare(Initial_Sample_Size,"2,350 individuals in 549 families; 390 trios")=0 then Initial_Sample_Size="2,350 individuals in 549 families";
if compare(Initial_Sample_Size,"200 > 85th pct,200 < 15th pct,7,817 cohort members")=0 then Initial_Sample_Size="200; 200; 7,817";
if compare(Initial_Sample_Size,"5,312-9,707 individuals")=0 then Initial_Sample_Size="9707";
if compare(Initial_Sample_Size,"7,440-10,783 individuals")=0 then Initial_Sample_Size="10783";


*Exclusions;

*1.Race;
Initial_Sample_Size=lowcase(Initial_Sample_Size);
Initial_Sample_Size=compress(Initial_Sample_Size,",");

Replication_Sample_Size=lowcase(Replication_Sample_Size);
Replication_Sample_Size=compress(Replication_Sample_Size,",");

if find(Initial_Sample_Size,"japanese")>0 then race_del=1;
else if find(Initial_Sample_Size,"chinese")>0 then race_del=1;
else if find(Initial_Sample_Size,"african")>0 then race_del=1;
else if find(Initial_Sample_Size,"filipino")>0 then race_del=1;
else if find(Initial_Sample_Size,"micronesian")>0 then race_del=1;
else if find(Initial_Sample_Size,"indian")>0 then race_del=1;
else if find(Initial_Sample_Size,"asian")>0 then race_del=1;
else if find(Initial_Sample_Size,"korean")>0 then race_del=1;
else if find(Initial_Sample_Size,"pima indians")>0 then race_del=1;
else if find(Initial_Sample_Size,"kosraen")>0 then race_del=1;

* BUT CAUCASIAN IS GOOD, so Keep it;
if find(Initial_Sample_Size,"caucasian")>0 then race_del=.;
if find(Initial_Sample_Size,"european")>0 then race_del=.;

* THIS PAPER IS GOOD AND WE KEEP it;
if First_Author="Elliott P" and Journal="JAMA" then race_del=.;
if First_Author="Kilpeläinen TO" and Journal="Nat Genet" then race_del=.;
if First_Author="Smith NL" and Journal="Circ Cardiovasc Genet" then race_del=.;

* THIS P_value is wrong and we need to keep it;
*if First_Author="Schunkert H" and SNPs="rs3184504" then p_Value="6E-6";

* Delete not european races;
if race_del=1 then delete;

*2. P-value;
pval=input(p_Value,best12.);
if pval > 5E-8 or  pval=.  /*or RA in ("?","N","P")*/ then delete;

run;


proc freq data=step2 noprint;
table Initial_Sample_Size  / out=temp;
run;

*3. Sample size;
    * Subjects with initial sample size + replication sample size <10000 are dropped;
data step3;
set step2;


array start(5) start1-start5;
array stop(5) stop1-stop5;
array initial(5) initial1-initial5;

array startr(5) startr1-startr5;
array stopr(5) stopr1-stopr5;
array repl(5) repl1-repl5;

if Replication_Sample_Size="nr" then Replication_Sample_Size=0;

do i=1 to 5;
*Initial;
if i=1 then start[i]=anydigit(Initial_Sample_Size,start[i]);
else start[i]=anydigit(Initial_Sample_Size,stop[i-1]+1);
stop[i]=notdigit(Initial_Sample_Size,start[i]);
initial[i]=input(substrn(Initial_Sample_Size,start[i],stop[i]-start[i]),best12.);
*some conditions;
if initial[i]=. then initial[i]=0;
if i gt 1 and initial[i-1]=0 then initial[i]=0;

*Replication;
if i=1 then startr[i]=anydigit(Replication_Sample_Size,startr[i]);
else startr[i]=anydigit(Replication_Sample_Size,stopr[i-1]+1);
stopr[i]=notdigit(Replication_Sample_Size,startr[i]);
repl[i]=input(substrn(Replication_Sample_Size,startr[i],stopr[i]-startr[i]),best12.);
*some conditions;
if repl[i]=. then repl[i]=0;
if i gt 1 and repl[i-1]=0 then repl[i]=0;

end;

initi=initial1+initial2+initial3+initial4+initial5;
rep=repl1+repl2+repl3+repl4+repl5;

*Exclusion;
if initi+rep<10000 then delete;

RA_freq = input(Risk_Allele_Frequency,best12.);
RA=input(compress(Strongest_SNP_Risk_Allele,"1234567890rs-"),$1.);
if RA='?' then RA=.;

*Databse weights;
w=input(OR_or_beta,best12.);

*Direction;
direction="";

drop i start1-start5 stop1-stop5 startr1-startr5 stopr1-stopr5 initial1-initial5 repl1-repl5 Journal Link 
      Platform__SNPs_passing_QC_ CNV Risk_Allele_Frequency p_Value OR_or_beta Strongest_SNP_Risk_Allele race_del
      Chr_pos Region Upstream_gene_id Downstream_gene_id Snp_gene_ids Upstream_gene_distance Downstream_gene_distance
	  Pvalue_mlog Context Intergenic; 


run;

* Count number of SNPs;
proc freq data=step3 noprint;
table SNPs/ out=temp;
run;

proc freq data=step3 noprint;
table Study/ out=temp;
run;

proc freq data=step3 noprint;
table trait/ out=temp;
run;

** CHECK ANOTHER TIME RACE ;
** CHECK ANOTHER TIME SAMPLE SIZES;

proc freq data=step3;
table Initial_Sample_Size Replication_Sample_Size ;
run;


*****************************
3. Add CHD for sensitivity  *
*****************************;

* Add teslovich (see later why);
PROC IMPORT OUT= WORK.lipids 
            DATAFILE= "&path.\data\SNPs_selection\Lipids.txt" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data lipids; set lipids; if SNPs="rs495828" then delete; NHGRI=1; run;


***************************************************
4. MERGE, CREATE GROUPS AND SOLVE MINOR PROBLEMS  *
***************************************************;


*** First print all the p_Value__text_ to see if we need to rename traits;
proc freq data=step3;
table p_Value__text_ / out=temp; 
run;

*** Check document: /Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/SNPs_selection/Checks/P_value_text.xlxs;

proc format;

value grFHSf 1= "TC"
2="HDL-C"
3="Smoke"
4="SBP"
5="BMI"
6="diabetes"
7="CHD";

run;


data step4;
set step3 lipids;

*Fix some SNPs problem;
if compare(SNPs,"rs2048327,rs3")=0 then SNPs="rs2048327";
if compare(SNPs,"rs1531343a")=0 then SNPs="rs1531343";


*Recode traits to better trait;
if find(p_Value__text_,"(Waist)")>0 then trait=78;
if find(p_Value__text_,"(WHR in women)")>0 then trait=80;
if find(p_Value__text_,"(WC)")>0 then trait=78;
if find(p_Value__text_,"(TRIG)")>0 then trait=70;
if find(p_Value__text_,"(UA)")>0 then trait=115;
if find(p_Value__text_,"(TG)")>0 then trait=70;
if find(p_Value__text_,"(SBP)")>0 then trait=69;
if find(p_Value__text_,"(QT interval)")>0 then trait=60;
if find(p_Value__text_,"(PR interval)")>0 then trait=58;
if find(p_Value__text_,"(LDL)")>0 then trait=46;
if find(p_Value__text_,"(HR)")>0 then trait=62;
if find(p_Value__text_,"(HDL)")>0 then trait=35;
if find(p_Value__text_,"(GLU)")>0 then trait=32;
if find(p_Value__text_,"(FPG)")>0 then trait=32;
if find(p_Value__text_,"(CRP)")>0 then trait=15;
if find(p_Value__text_,"(AngCAD)")>0 then trait=21;
if find(p_Value__text_,"(AngCAD/MI)")>0 then trait=21;


if compare(p_Value__text_,"(Gout)")=0 then delete;
if compare(p_Value__text_,"(CCTC)")=0 then delete;
if compare(p_Value__text_,"(CTTG)")=0 then delete;



if compare(p_Value__text_,"(FPG)")=0 then _5__CI__text_="increase";
if compare(p_Value__text_,"(FI)")=0 then _5__CI__text_="increase";


* Now delete the outcomes not needed anymore;
if p_Value__text_ in ("(Waist)","(WHR in women)","(WC)","(TRIG)","(UA)","(TG)","(SBP)","(QT interval)","(PR interval)","(LDL)",
"(HR)","(HDL)","(GLU)","(FPG)","(CRP)","(AngCAD)","(AngCAD/MI)","(women)","(whites)","(obese)","(non-obese)","(men)","(ischemic stroke)",
"(children)","(DGI+FUSION+WTCCC)","(EA)","(Urate)","(discovery + validation)") then p_Value__text_=.;


*More stricht groups traits (FHS + BMI);
if trait in (19) then grFHS = 1;
else if  trait in (35) then grFHS = 2;
else if trait in (63) then grFHS = 3;
else if trait in (69) then grFHS = 4;
else if trait in (12,13,14,54,55,56,122,123) then grFHS = 5;
else if trait in (23,73,74,75) then grFHS = 6;
else if trait in (21,52,53) then grFHS = 7;

if NHGRI=. then NHGRI=0;


format grFHS grFHSf.;

drop Initial_Sample_Size Replication_Sample_Size;
run;


proc freq data=step4 noprint;
table trait/ out=temp;
run;

****************************************************
5. CHECK: MISSING RISK ALLELE OR ALLELE FREQUENCY  *
****************************************************;

* View file: 

* Check missing risk allele or missing allele frequency;
data problem;
set step5;
if RA="" or RA="." or RA_freq=.;
run;


** MANUALLY ADD MISSING Risk Allele or Allele frequency;

data step5;
set step4;


* Barrett needs to be excluded. We couldn't understand which is the minor allele (to associate the minor allele frequency);
if First_Author="Barrett JC" then delete;

else if RA="." and First_Author="Bradfield JP" and SNPs="rs478222" then RA="A";
else if RA="." and First_Author="Bradfield JP" and SNPs="rs539514" then RA="T";
else if RA="." and First_Author="Bradfield JP" and SNPs="rs924043" then RA="C";


else if RA="." and First_Author="Bradfield JP" and SNPs="rs9299" then do;
RA="T";
RA_freq=0.652;
end;
else if RA="." and First_Author="Bradfield JP" and SNPs="rs9568856" then do;
RA="A";
RA_freq=0.158;
end;

else if RA_freq=. and First_Author="Chambers JC" then delete;

else if RA_freq=. and First_Author="Cooper JD" then delete;

else if RA_freq=. and First_Author="Dastani Z" and SNPs="rs12051272" then RA_freq=0.03;
else if RA_freq=. and First_Author="Dastani Z" and SNPs="rs182052" then delete;
else if RA_freq=. and First_Author="Dastani Z" and SNPs="rs6488898" then delete;

else if RA_freq=. and First_Author="Davies RW" and SNPs="rs3869109" then RA_freq=0.55;

else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs2794520" then RA_freq=0.66;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs12037222" then RA_freq=0.24;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs1183910" then RA_freq=0.67;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs4420065" then RA_freq=0.61;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs4129267" then RA_freq=0.60;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs1260326" then RA_freq=0.41;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs6734238" then RA_freq=0.42;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs9987289" then RA_freq=0.90;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs340029" then RA_freq=0.62;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs13233571" then RA_freq=0.86;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs12239046" then RA_freq=0.61;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs10745954" then RA_freq=0.50;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs1800961" then RA_freq=0.95;
else if RA_freq=. and First_Author="Dehghan A" and SNPs="rs10521222" then RA_freq=0.94;

else if RA="." and First_Author="Dehghan A" and SNPs="rs1165205" then RA="A";
else if RA="." and First_Author="Dehghan A" and SNPs="rs2231142" then RA="T";

else if RA_freq=. and First_Author="Elliott P" and SNPs="rs6700896" then RA_freq=0.38;
else if RA_freq=. and First_Author="Elliott P" and SNPs="rs1183910" then RA_freq=0.32;
else if RA_freq=. and First_Author="Elliott P" and SNPs="rs4537545" then RA_freq=0.43;
else if RA_freq=. and First_Author="Elliott P" and SNPs="rs4420638" then RA_freq=0.19;
else if RA_freq=. and First_Author="Elliott P" and SNPs="rs7553007" then RA_freq=0.33;

else if RA_freq=. and First_Author="Grant SF" then RA_freq=0.474;

else if RA="." and First_Author="Gudbjartsson DF" then delete;
else if RA_freq=. and First_Author="Gudbjartsson DF" then RA_freq=0.183;

else if RA_freq=. and First_Author="Heard-Costa NL" and SNPs="rs10146997" then RA_freq=0.21;
else if RA_freq=. and First_Author="Heard-Costa NL" and SNPs="rs1558902" then delete;

else if RA="." and First_Author="Huang J" and SNPs="rs2227631" then delete;
else if RA="." and First_Author="Huang J" and SNPs="rs6486122" then delete;
else if RA="." and First_Author="Huang J" and SNPs="rs6976053" then delete;


else if RA="." and First_Author="Huang J" and SNPs="rs1265564" then do;
RA="A";
RA_freq=0.48;
end;
else if RA="." and First_Author="Huang J" and SNPs="rs61839660" then do;
RA="C";
RA_freq=0.09;
end;
else if RA="." and First_Author="Huang J" and SNPs="rs7018475" then do;
RA="T";
RA_freq=0.24;
end;

else if First_Author="Kraja AT" then delete;

else if RA="." and First_Author="Loos RJ" then delete;

else if First_Author="Mehta NN" then delete;

else if First_Author="Mitchell GF" then delete;

else if RA_freq=. and First_Author="O'Donnell CJ" and SNPs="rs9349379" then RA_freq=0.59;
else if RA_freq=. and First_Author="O'Donnell CJ" and SNPs="rs1333049" then RA_freq=0.47;

else if RA_freq=. and First_Author="Paternoster L" then RA_freq=0.17;

else if RA="." and First_Author="Qi L" then do;
  RA_freq=0.23;
  RA="T";
end;

else if RA_freq=. and First_Author="Saxena R" and SNPs="rs1260326" then RA_freq=0.40;
else if RA_freq=. and First_Author="Saxena R" and SNPs="rs2877716" then RA_freq=0.77;
else if RA_freq=. and First_Author="Saxena R" and SNPs="rs10423928" then RA_freq=0.18;

else if RA="." and First_Author="Smith NL" then RA="G";

else if RA="." and First_Author="Steinthorsdotti" then delete;

else if RA="." and First_Author="Surakka I" then delete;

else if First_Author="The Coronary Ar" and SNPs="rs9349379" then delete;
else if First_Author="The Coronary Ar" and SNPs="rs646776" then delete;
else if First_Author="The Coronary Ar" and SNPs="rs4977574" then delete;

else if RA="" and First_Author="Teslovich TM" and SNPs="rs174546" then RA="T";

else if RA="." and First_Author="Timpson NJ" then delete;


else if RA_freq=. and First_Author="Voight BF" and SNPs="rs1387153" then RA_freq=0.28;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs7578326" then RA_freq=0.64;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs1531343" then RA_freq=0.10;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs13292136" then RA_freq=0.93;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs243021" then RA_freq=0.46;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs8042680" then RA_freq=0.22;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs11634397" then RA_freq=0.60;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs7957197" then RA_freq=0.85;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs1552224" then RA_freq=0.88;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs231362" then RA_freq=0.52;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs896854" then RA_freq=0.48;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs4457053" then RA_freq=0.26;
else if RA_freq=. and First_Author="Voight BF" and SNPs="rs972283" then RA_freq=0.55;
else if RA_freq=. and First_Author="Voight BF" then delete; *delete the other SNPs;	


else if RA_freq=. and First_Author="Zeggini E" and SNPs="rs6931514" then delete;

run;


****************************************************************************
6. ADD PROXIES FOR SNPs not FOUND in HAPMAP relase 27 and in TWINGENE GWAS *
****************************************************************************;

** Some SNPs cannot be found in HapMap relase 27 and I need to use proxies;
** To do that I use http://www.broadinstitute.org/mpg/snap/ldsearch.php;
/* SNPs LIST:
rs1558861
rs2779116
rs3008621
rs3812316
rs4836133
rs9411489
*/
* Check also the proxy-xls documents;

** One SNPs have been merged and i used the update rsID, but this is not available so deleted;
*rs7826222;

data step6;
set step5;

** THIS SNPs are not in HapMap release 27 and I need to use proxyes;
if SNPs="rs1558861" then do;
	SNPs="rs6589566"; RA="G"; RA_freq=0.07; proxy=1;
end;

if SNPs="rs2779116" then do;
	SNPs="rs2518491"; RA="T"; RA_freq=0.32;proxy=1;
end;

if SNPs="rs3008621" then do;
	SNPs="rs4240934"; RA="C"; RA_freq=0.90;proxy=1;
end;

if SNPs="rs3812316" then do;
	SNPs="rs17145738"; RA="T"; RA_freq=0.12;proxy=1;
end;

if SNPs="rs4836133" then do;
	SNPs="rs6864049"; RA="G"; RA_freq=0.49;proxy=1;
end;

if SNPs="rs9411489" then do;
	SNPs="rs495828"; RA="T"; RA_freq=0.19;proxy=1;
end;

if SNPs="rs326" then do;
	SNPs="rs13702"; RA="T"; RA_freq=0.712;proxy=1;
end;

if SNPs="rs28927680" then do;
	SNPs="rs17120029"; RA="T"; RA_freq=0.04;proxy=1;
end;

* I can not find proxyes;
if SNPs="rs7826222" then delete; * merged and not available;
if SNPs="rs3798220" then delete; * This in not included in the CHIP and not proxy are available;
if SNPs="rs1265564" then delete;	*1000genome imputation;
if SNPs="rs61839660" then delete;	*1000genome imputation;

run;






***************
7. ANNOTATION *
***************;


** EXPORT SINGLE SNPS TO BE ANNOTATED in http://hapmap.ncbi.nlm.nih.gov/biomart/martview/3e9af05f2a67fe20a2302eb8c3a8a5b3;
* flag these options:  chromosome,position,marker id,alleles,reference allele,reference allele frequency,other allele,other allele frequency;
proc sort data=step6 nodupkey out=exp (keep=SNPs); by SNPs; run;

PROC EXPORT DATA= WORK.exp
            OUTFILE= "&path.\Data\SNPs_selection\allSNPs_tobe_annotated.txt" 
            DBMS=TAB REPLACE;
     PUTNAMES=NO;
RUN;



*####################################*
*####################################*
* Pgm breaks, manual action required *
*####################################*
*####################################*


** IMPORT ANNOTATED snps;
PROC IMPORT OUT= WORK.EXP_ANNOTATED 
            DATAFILE= "&path.\Data\SNPs_selection\allSNPs_annotated.txt" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;


data exp_annotated; set exp_annotated; rename marker_id=SNPs; run;
proc sort data=exp_annotated; by SNPs; run;
proc sort data=step6; by SNPs; run;




*******************************
8. CHECK: ANNOTATION PROBLEMS *
*******************************;
* Chck this folder for more info: /Users/AndreaGanna/Documents/Work/Phd_KI/Genscore/Documents/SNPs_selection/Checks/Annotation_problems;


data def;
merge exp_annotated step6;
by SNPs;


if alleles in  ("C/T","A/G","A/C","G/T") then do;

     ** 1. if our allele correspond to one of the HapMap alleles then should be on + strand 
           (check anyway for too big differences in allele frequencies);

    if RA=reference_allele or RA=other_allele then check=1;

     ** 2. if our allele correspond to one of the flipped HapMap alleles then should be on - strand 
           (check anyway for too big differences in allele frequencies);
	else if (RA="T" and reference_allele="A" or RA="A" and reference_allele="T") 
         or (RA="T" and other_allele="A" or RA="A" and other_allele="T") 
		 or (RA="C" and reference_allele="G" or RA="G" and reference_allele="C")
		 or (RA="C" and other_allele="G" or RA="G" and other_allele="C") then check=2;	

end;


*** DIFFERENCES;
if check=1 then do;
  if RA=reference_allele then diff=abs(reference_allele_frequency-RA_freq);
  if RA=other_allele then diff=abs(other_allele_frequency-RA_freq);
end;


if check=2 then do;
	     if  RA="T" and reference_allele="A" then diff = abs(reference_allele_frequency-RA_freq);
         if  RA="A" and reference_allele="T" then diff = abs(reference_allele_frequency-RA_freq);
         if  RA="T" and other_allele="A" then diff = abs(other_allele_frequency-RA_freq);
         if  RA="A" and other_allele="T" then diff = abs(other_allele_frequency-RA_freq);
		 if  RA="C" and reference_allele="G" then diff = abs(reference_allele_frequency-RA_freq);
         if  RA="G" and reference_allele="C" then diff = abs(reference_allele_frequency-RA_freq);
		 if  RA="C" and other_allele="G" then diff = abs(other_allele_frequency-RA_freq);
         if  RA="G" and other_allele="C" then diff = abs(other_allele_frequency-RA_freq);
end;





if alleles in  ("A/T","C/G") then do;

   ** 3. Alleles A/T and C/G with allele frequencies around 0.50 are instable and need to be manually checked;
   if 0.40<RA_freq<0.60 then do;
     if RA=reference_allele or RA=other_allele then check=3;
   end;

   else do;

   ** 4. if our allele correspond to one of the HapMap alleles and the allele frequency is almost the same 
          then should be on + strand;
   if (RA=reference_allele  and RA_freq-0.10 < reference_allele_frequency < RA_freq+0.10)
      or (RA=other_allele and RA_freq-0.10 < other_allele_frequency < RA_freq+0.10) then check=4;

   ** 5. On - strand or any other kind of error;
   else  check=5;
   end;

end;

*** DIFFERENCES;
if check=4 then do;
  if RA=reference_allele then diff=abs(reference_allele_frequency-RA_freq);
  if RA=other_allele then diff=abs(other_allele_frequency-RA_freq);
end;


run;


* Check problem with allele alligment with HapMap;
** That means to check big differences and check=5,3;
data problem;
set def;
if diff ge 0.1 or check in (5,3);
run;


data def;
set def;

** Adjusting problems (check document check.xls for more info (is in the last_quality folder) );
if first_author="Heid IM" then delete;  *** THE STUDY HAS SOME PROBLEMS THAT I COULDN'T SOLVE ****;
else if first_author="Cooper JD" and SNPs="rs3825932" then check=1;
else if first_author="Holm H" and SNPs="rs1321311" then check=2;
else if first_author="Huang J" and SNPs="rs7018475" then delete;
else if first_author="Kathiresan S" and SNPs="rs4977574" then check=1;
else if first_author="Kooner JS" and SNPs="rs2075292" then check=1;
else if first_author="Meyre D" and SNPs="rs1424233" then check=2;
else if first_author="Saxena R" and SNPs="rs13266634" then check=1;
else if first_author="Thorleifsson G" and SNPs="rs6499640" then check=1;
else if first_author="Waterworth DM" and SNPs="rs442177" then check=2;
else if first_author="Waterworth DM" and SNPs="rs174548" then check=1;
else if first_author="Erdmann J" and SNPs="rs3739998" then check=2;
else if first_author="Richards JB" and SNPs="rs1648707" then delete;
else if first_author="Richards JB" and SNPs="rs4311394" then delete;
else if First_Author="Dehghan A" and SNPs="rs9987289" then RA="G";
else if First_Author="Dupuis J" and SNPs="rs11558471" then RA_freq=0.68;
else if first_author="The Coronary Ar" and SNPs="rs974819" then check=1;

** All SNPs with check = 3 need to be considered as check =1;

* Teslovich is removed and added manualy;
else if first_author="Teslovich TM" and Snp_id_current ne . then delete;
if SNPs="rs9411489" then delete;

** Flip strand (check=2);
if check=2 then do;
  if RA="G" then RA="C";
  else if RA="C" then RA="G";
  else if RA="A" then RA="T";
  else if RA="T" then RA="A";
end;
run;



** Check the traits;

proc freq data=def noprint;
table trait/ out=temp;
run;


******************************
9. FLIP - DEFINITIVE DATASET *
******************************;

data def_sin;
set def;

**1) Exclude traits where we are not sure of the direction or have a U shape;
if trait in (95,60) or p_Value__text_ in ("(QRS complex)","(LV internal diastolic dime","(HOMA-B)") then delete;


** 2) Assign direction to dichotomous traits (all risk factors);
if  trait in (1,7,8,21,36,39,45,52,53,54,56,64,68,72,73,75,117) then do;
	dicotomous=1;
	direction="+";
end;
* Hypertension;
if First_Author="Ehret GB" and trait=39 then do;
	if find(_5__CI__text_,"increase")>0 then direction="+";
	else if find(_5__CI__text_,"decrease")>0 then direction="-";
end;

* Type 1 diabetes;
if First_Author="Bradfield JP" and trait=72 then direction="+";


** 3) Assign direction to non-dichotomous traits;
if  trait not in (1,7,8,21,36,39,45,52,53,54,56,64,68,72,73,75,117)	then do;
	if find(_5__CI__text_,"increase")>0 then direction="+";
	else if find(_5__CI__text_,"decrease")>0 then direction="-";
	else if find(_5__CI__text_,"higher")>0 then direction="+";
	else if find(_5__CI__text_,"lower")>0 then direction="-";
end;

* Obesity (extreme);
if First_Author="Paternoster L" and SNPs="rs734597" then direction="+";

* Remove this SNP because effec tonly in women;
if First_Author="Fox CS" and SNPs="rs1659258" then delete;
if First_Author="Fox CS" and SNPs="rs11118316" then direction="+";



* (HOMA-IR);
if p_Value__text_ in ("(HOMA-IR)") then direction="+";



if First_Author="Thorleifsson G" then direction="+";
if First_Author="Willer CJ" and SNPs="rs10938397" then direction="+"; 
if First_Author="Frayling TM" and SNPs="rs9939609" then direction="+"; 

* Proinsuline levels;
if First_Author="Strawbridge RJ" then direction="+";


** 4) extra problematic directions;
if first_author="The Tobacco and" and SNPs="rs3025343" then direction="-"; *smoking cessation;

*Error in the catalog, the RA is T, moreover T decrease the risk of smoke initiation then is protective, so direction is -;
else if first_author="The Tobacco and" and SNPs="rs6265" then do;
  direction="-"; 
  RA="T";
end;

** 4) HDL coholesterol and adiponectin have inverted directions;
if trait in (35,2) and direction="+" then direction="-";
else if trait in (35,2) and direction="-" then direction="+";



*** FINALIZE;
**1) Determine if risk allele = to hapmap reference allele or not;
if RA=reference_allele then rev=1;
else if RA=other_allele then rev=0;

*** 2) And now revert allele if direction = -;
if direction="-" then do;
	if rev=1 then RA=other_allele;
	else if rev=0 then RA=reference_allele;
	RA_freq=1-RA_freq;
end;

w=abs(w);


if SNPs="rs495828" then delete;
run;




********************************************
10. DUPLICATED SNPs within/between  TRAITS *
********************************************;

******** WITHIN ********;

proc sort data=def_sin; by trait; run;
proc freq data=def_sin noprint;
table SNPs / out=temp;
by trait;
run;



proc sort data=def_sin; by trait SNPs; run;
proc sort data=temp; by trait SNPs; run;
data dupl_trait_in;
merge temp def_sin;
by trait SNPs;
keep SNPs trait RA COUNT First_Author Study RA_freq alleles;
if COUNT > 1;
run;



* Keep only duplicated with non-concordant RA;
proc sort data=dupl_trait_in; by index SNPs; run;
data dupl_trait_in_f;
retain RA2;
set dupl_trait_in; by trait SNPs;
if first.SNPs then RA2=RA;
if RA2=RA then c=1;
if c= .;
drop RA RA2 c trait First_Author Study RA_freq alleles;
run;

proc sort data=dupl_trait_in_f; by SNPs; run;
proc sort data=dupl_trait_in; by SNPs; run;
data dupl_trait_in_f2;
merge dupl_trait_in dupl_trait_in_f (in=x);
by SNPs;
if x=1;
run;


** Adjust Def;
data def_sin;
set def_sin;
if first_author="Kooner JS" and SNPs="rs17145738" then delete; 
if first_author="Kathiresan S" and SNPs="rs6511720" then delete;
run;


******** BETWEEEN ********;



proc sort data=def_sin out=no_dup; by trait SNPs descending initi;
run;


data no_dup;
set no_dup;
by trait SNPs descending initi;
if first.SNPs; 
run;

*** TO IDENTIFY DUPLICATED SNPs WITHIN DIFFERENTS TRAITS WHEN THEY HAVE DISCORDANT RA;

proc freq data=no_dup noprint;
table SNPs / out=temp;
run;

proc sort data=no_dup; by SNPs; run;
proc sort data=temp; by SNPs; run;
data dupl_trait_bt;
merge temp no_dup;
by SNPs;
keep SNPs trait RA COUNT First_Author Study RA_freq alleles;
if COUNT > 1;
run;

* Keep only duplicated with non-concordant RA;
proc sort data=dupl_trait_bt; by SNPs; run;
data dupl_trait_bt_f;
retain RA2;
set dupl_trait_bt; by SNPs;
if first.SNPs then RA2=RA;
if RA2=RA then c=1;
if c= .;
drop RA RA2 c trait First_Author Study RA_freq alleles;
run;

data dupl_trait_bt_f2;
merge dupl_trait_bt dupl_trait_bt_f (in=x);
by SNPs;
if x=1;
run;



** Adjust Def_sin;


data def_sin;
set def_sin;

** Keep only the RA that we want;
if first_author="Speliotes EK" and SNPs="rs13107325" then delete;
else if first_author="Teslovich TM" and SNPs="rs13107325" then delete;
else if first_author="Eijgelsheim M" and SNPs="rs174547" then delete;
else if first_author="Middelberg RP" and SNPs="rs1800562" then delete;
else if first_author="Dehghan A" and SNPs="rs1800961" then delete;
else if first_author="Middelberg R" and SNPs="rs2075650" then delete;
else if first_author="Dehghan A" and SNPs="rs4420638" then delete;
else if first_author="Dupuis J" and SNPs="rs4607517" and RA="G" then delete;
else if first_author="Willer CJ" and SNPs="rs4775041" then delete;
else if first_author="Aulchenko YS" and SNPs="rs4939883" and RA="C" then delete;
else if first_author="Kathiresan S" and SNPs="rs6511720" then delete;
else if first_author="Aulchenko YS" and SNPs="rs780094" then delete;
else if first_author="Kolz M" and SNPs="rs780094" then delete;
else if first_author="Dehghan A" and SNPs="rs9987289" then delete;
else if first_author="Pfeufer A" and SNPs="rs3807989" then delete;


run;



*************************
11. SNPs quality metrics*
*************************;
PROC IMPORT OUT= WORK.twge_inf 
            DATAFILE= "&path.\Data\inf\twge_inf" 
            DBMS=DLM REPLACE;
     DELIMITER='20'x; 
     GETNAMES=YES;
     DATAROW=2; 
RUN;
data twge_inf; set twge_inf; SNPs=rs_id; info_twge=info; keep info_twge SNPs; run;


PROC IMPORT OUT= WORK.ulsam_inf 
            DATAFILE= "&path.\Data\inf\ulsam_inf" 
            DBMS=DLM REPLACE;
     DELIMITER='20'x; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data ulsam_inf; set ulsam_inf; SNPs=VAR2; info_ulsam=VAR5; keep info_ulsam SNPs; run;
/*
PROC IMPORT OUT= WORK.ulsam_inf2 
            DATAFILE= "&path.\Data\inf\ulsam_inf2.txt" 
            DBMS=DLM REPLACE;
	     DELIMITER='20'x; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data ulsam_inf2; set ulsam_inf2; SNPs=VAR2; info_ulsam=VAR5; keep info_ulsam SNPs; run;
*/

PROC IMPORT OUT= WORK.gosh_inf 
            DATAFILE= "&path.\Data\inf\gosh_inf" 
            DBMS=DLM REPLACE;
     DELIMITER='20'x; 
     GETNAMES=NO;
     DATAROW=1; 
RUN;
data gosh_inf; set gosh_inf; SNPs=VAR2; info_gosh=VAR5; keep info_gosh SNPs; run;


proc sort data=twge_inf; by SNPs;run;
proc sort data=ulsam_inf; by SNPs;run;
proc sort data=gosh_inf; by SNPs;run;
proc sort data=def_sin; by SNPs;run;

data temp;
merge def_sin (in=x) ulsam_inf twge_inf gosh_inf;
by SNPs;
if x=1;
run;

data def_sin; set temp;	
if info_ulsam < 0.4 and info_gosh < 0.4 and info_twge < 0.4 then delete;
run;

data low_imp;
set def_sin;
where info_ulsam < 0.4 or  info_gosh < 0.4; 
keep SNPs;
run;  

proc sort data=low_imp nodupkey out=low_imp_exp; by SNPs; run;

PROC EXPORT DATA= WORK.low_imp_exp
            OUTFILE= "&path.\Results\SNPs_selection\low_imp.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;



* Count number of SNPs;
proc sort data=def_sin nodupkey out=count; by SNPs; where NHGRI=1;run;
proc freq data=count noprint;
table SNPs/ out=temp;
where NHGRI=1;
run;

proc freq data=count noprint;
table Study/ out=temp;
where NHGRI=1;
run;

proc freq data=count noprint;
table trait/ out=temp; 
where NHGRI=1;
run; 

********************
12. EXPORT SNPLIST *
********************;

proc sort data=def_sin out=snplist; by trait SNPs descending initi;
run;


*** TO DROP DUPLICATED SNPs INSIDE SAME TRAIT (keeping the one with highest sample size);
data snplist;
set snplist;
by trait SNPs descending initi;
if first.SNPs; 
run;


proc sort data=snplist nodupkey out=snplist_e; by SNPs; run;
data snplist_e;
informat chromosome SNPs RA pval;
set snplist_e;
if chromosome="chrX" then delete;
keep  chromosome SNPs RA pval;
run;

PROC EXPORT DATA= WORK.snplist_e
            OUTFILE= "&path.\Results\SNPs_selection\snplist.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


*********************************************
e1. EXPORT: CATALOG SCORE (ALLCAT) NHGRI=1  *
*********************************************;
proc sort data=def_sin out=allcat; by trait SNPs descending initi; where NHGRI=1;
run;


*** TO DROP DUPLICATED SNPs INSIDE SAME TRAIT  (keeping the one with highest sample size);
data allcat;
set allcat;
by trait SNPs descending initi;
if first.SNPs; 
run;

proc sort data=allcat nodupkey out=allcat_e; by SNPs; run;
data allcat_e;
informat chromosome SNPs RA pval;
set allcat_e;
if chromosome="chrX" then delete;
keep  chromosome SNPs RA pval;
run;

PROC EXPORT DATA= WORK.allcat_e
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;



*********************************************
e2. EXPORT: FHS specific trait from catalog *
********************************************;

** BMI;
data exp_bmi;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=5;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_bmi nodupkey out=exp_bmi; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_bmi
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_bmi.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;

** SBP;
data exp_sbp;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=4;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_sbp nodupkey out=exp_sbp; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_sbp
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_sbp.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


** smoke;
data exp_smoke;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=3;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_smoke nodupkey out=exp_smoke; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_smoke
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_smoke.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


** HDL;
data exp_hdl;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=2;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_hdl nodupkey out=exp_hdl; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_hdl
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_hdl.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;

** TC;
data exp_tc;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=1;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_tc nodupkey out=exp_tc; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_tc
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_tc.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


** T2d;
data exp_diab;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=6;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_diab nodupkey out=exp_diab; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_diab
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_t2d.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


** chd;
data exp_chd;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS=7;
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_chd nodupkey out=exp_chd; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_chd
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_chd.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;


** Export all FHS trait;
data exp_fhs;
informat chromosome SNPs RA pval;
set def_sin;
keep chromosome SNPs RA pval;
if grFHS in (1,2,3,4,5,6);
if chromosome="chrX" then delete;
where NHGRI=1;
run;

proc sort data=exp_fhs nodupkey out=exp_fhs; by SNPs; run;

PROC EXPORT DATA= WORK.Exp_fhs
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_fhs.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;



***********************************************************
e5. EXPORT: Annotation information for supplemetary table *
***********************************************************;


proc sort data=def_sin out=allcat; by trait SNPs descending initi; where NHGRI=1;
run;


*** TO DROP DUPLICATED SNPs INSIDE SAME TRAIT  (keeping the one with highest sample size);
data allcat;
set allcat;
by trait SNPs descending initi;
if first.SNPs; 
if first_author="Teslovich TM" then PUBMEDID="20686565";
if first_author="Deloukas" then PUBMEDID="23202125";
if first_author="Morris" then PUBMEDID="22885922";
keep SNPs Reported_Gene_s_ PUBMEDID Disease_Trait;
run;

proc sort data=allcat nodupkey out=allcat_e; by SNPs; run;


PROC EXPORT DATA= WORK.allcat_e
            OUTFILE= "&path.\Results\SNPs_selection\annotation.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=YES;
RUN;



****************************************************************************
e_review. EXPORT: CATALOG SCORE (ALLCAT) NHGRI=1, but no CHD-related SNPs  *
****************************************************************************;
proc sort data=def_sin out=allcat; by trait SNPs descending initi; where NHGRI=1 and grFHS ne 7;
run;


*** TO DROP DUPLICATED SNPs INSIDE SAME TRAIT  (keeping the one with highest sample size);
data allcat;
set allcat;
by trait SNPs descending initi;
if first.SNPs; 
run;

proc sort data=allcat nodupkey out=allcat_e; by SNPs; run;
data allcat_e;
informat chromosome SNPs RA pval;
set allcat_e;
if chromosome="chrX" then delete;
keep  chromosome SNPs RA pval;
run;

PROC EXPORT DATA= WORK.allcat_e
            OUTFILE= "&path.\Results\SNPs_selection\ALLCAT_no_chd.txt" 
            DBMS=dlm REPLACE;
	delimiter=' ';
     PUTNAMES=NO;
RUN;

