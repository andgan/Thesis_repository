PROJECTS NAME: Genscore - Multi-locus genetic risk scores for coronary heart disease 


- Import_TWINGENE_11JUL2012.sas, Import_ULSAM_11JUL2012.sas, Import_GOSH_05DEC2011.sas: to create final datasets for each study.

- CVD_GOSH_11JUL2012.sas: to attach the incident event outcomes to GOSH.

- SNP_selection_FINAL.sas: program to read the GWAS catalog and extract the final list of SNPs. 

- Grep_dosages_TWINGENE, Grep_dosages_ULSAM, Grep_dosages_GOSH. Bash script to elect the right SNPs dosages from the imputed data from each study.

- Flip_dosages_TWINGENE, Flip_dosages_ULSAM, Flip_dosages_GOSH. Bash script to flip the SNPs on the right risk allele. 

- Pruning_15JAN2013_review: bash and plink script to prune SNPs for each score.

- polygene_weights_15JAN2013: R script to create the polygene score

- Merge_flipped_dosages_and_score_FINAL_review: R script to merge all the scores from each study.

- merge_ULSAM_TWINGENE_GOSH_score_polygene_27SEP2011_review: put together the phenotypes from all the studies and also the genetic scores.
 
- ANA_15JAN2013_review: R script that can be used to analyze the results.