#----------------------------------------------
# Filename: extract_snps.py
# Study: metabo - CHD
# Author: Andrea Ganna
# Date: 27FEB2014
# Updated: 
# Purpose: Extract SNPs for CHD and candidate genes in Ulsam, Twingene and Pivus
# Note: 
#-----------------------------------------------
# Data used: twge.+.qced.b36.sample_1 impute2d.twge.+.qced.b36.chr$01$.ceu.p2.r22.int5M.controls.imp.gz snplist.txt
#            ULSAMb36.sample ULSAMb36_$01$.imputed.gz
#            pivus-omni-metabo.949.b36.sample PIVUSb36_$01$.imputed.gz
# Data created: exp_metabo_twge.csv exp_metabo_pivus.csv exp_metabo_ulsam.csv
#-----------------------------------------------
# OP: R 2.12.1
#-----------------------------------------------*/


## To run in uppmax
source .venvburrito/startup.sh

## To work on virtual enviroment
workon test


import pandas				
from pandas import *
import csv
import os
import sys

class get_snp_from_imp(object):
	
	def __init__(self, sample_file, imp_file, snplist_file):
		self.sample_file = sample_file
		self.imp_file = imp_file
		self.snplist_file = snplist_file
		
	def read_snps(self):
		with open(self.snplist_file) as snplist:
			self.snps = []
			self.chrs = []
			self.ras  = []
			num_lines = 0
			for line in snplist:
				num_lines += 1
				parts = line.split()
				if len(parts) >= 1 and len(parts) <=3 :
					self.chrs.append(parts[0])
					self.snps.append(parts[1])
					self.ras.append(parts[2])
			print "N. of SNPs read from snplist file: " + str(num_lines)
		   			
	def read_subj(self):
		with open(self.sample_file) as sample:
			next(sample)
			next(sample)
			self.ids = []
			num_lines = 0
			for line in sample:
				num_lines += 1
				parts = line.split()
				if len(parts) >= 2:
					self.ids.append(parts[0])
			print "N. of subjects read from the sample file: " + str(num_lines)		   
		   					
	def sum_dosages(self):
		dec = Series(self.parts[5:len(self.parts)]).apply(np.float32)
		res = Series([], dtype=float)
		i = 0
		while i < len(dec):
			a = Series(dec[i]*0+dec[(i+1)]*1+dec[(i+2)]*2, index=[i])
			res = res.append(a)
			i = i+3
		if self.r_all.upper() != self.parts[4] and self.r_all.upper() != self.parts[3]:
			print "Error with the risk allele"
		if self.r_all.upper() == self.parts[4]:
			self.res_n = res
			print "since RA = " + self.r_all + " and minor allele is " + self.parts[4] +": no flip"
		else:
			self.res_n = 2-res
			print "since RA = " + self.r_all + " and minor allele is " + self.parts[4] +": flip"
		return self.res_n	
			
	def read_imp(self):
		self.df = DataFrame(index=self.ids)
		print "Unzipping file"
		os.system("gunzip " + self.imp_file_m + " -c > .Temp")
		print "Unzipping finished"
		with 	open(".Temp") as imp:
			print "read chr" + str(chr)
			for line in imp:
				self.parts = line.split()
				if self.parts[1] in self.snps:
					print "read snp " + str(self.parts[1])
					# Find risk allele
					for r in range(len(self.ras)):
						if self.snps[r] == self.parts[1]:
							self.r_all = self.ras[r]					
					# Calculated correct dosages		
					self.sum_dosages()	
					self.res_n.index = self.ids
					self.df[self.parts[1] + "_" + self.r_all] = self.res_n
		os.system("rm .Temp")		
		return self.df			
	 
	def read_imp_loop(self):
		self.final_df = DataFrame(index=self.ids)
		for c in range(1,23):
			print "Read chromosome: " + str(c)
			if "chr"+str(c) in self.chrs:
				self.imp_file_m = self.imp_file.replace("$01$",str(c).zfill(2))
				self.read_imp()
				self.final_df = self.final_df.join(self.df)
			else:
				print "No SNPs in this chromosome"	
		return self.final_df
		
	def check_snps(self):
		snp_new=[]
		for e in range(len(self.snps)):
			snp_new.append(str(self.snps[e]) + '_' + self.ras[e])
		if set(self.final_df.columns) == set(snp_new):
			print 'All SNPs have been found'
		else:
			print str(list(set(snp_new) - set(self.final_df.columns))) + str(list(set(self.final_df.columns)-set(snp_new))) + 'are/is different'
				
				
def run_export(sample_file,imp_file,snplist_file,output_file):
	
	get_snp = get_snp_from_imp(sample_file,imp_file,snplist_file)
	get_snp.read_snps()
	get_snp.read_subj()
	out = get_snp.read_imp_loop()
	get_snp.check_snps()
	out['ID'] = out.index
	with open(output_file, 'wb') as csvfile:
		wr = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		wr.writerow(list(out.columns))
		for i in range(len(out.count(axis=1))):
			wr.writerow(out.ix[i,:])



run_export("/proj/b2011036/twge.gwas/twge.imp.hapmap.b36/impute.dosages/ceu.p2.r22/twge.+.qced.b36.sample_1","/proj/b2011036/twge.gwas/twge.imp.hapmap.b36/impute.dosages/ceu.p2.r22/impute2d.twge.+.qced.b36.chr$01$.ceu.p2.r22.int5M.controls.imp.gz" ,"/home/andrea/snplist.txt","/home/andrea/exp_metabo_twge.csv")

run_export("/proj/b2011036/ulsam.gwas/ulsam.imp.hapmap.b36/ULSAMb36.sample","/proj/b2011036/ulsam.gwas/ulsam.imp.hapmap.b36/ULSAMb36_$01$.imputed.gz" ,"/home/andrea/snplist.txt","/home/andrea/exp_metabo_ulsam.csv")

run_export("/proj/b2011036/pivus.gwas/pivus.imp.hapmap.b36/pivus-omni-metabo.949.b36.sample","/proj/b2011036/pivus.gwas/pivus.imp.hapmap.b36/PIVUSb36_$01$.imputed.gz" ,"/home/andrea/snplist.txt","/home/andrea/exp_metabo_pivus.csv")

