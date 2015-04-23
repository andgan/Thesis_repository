#! /usr/bin/Rscript

## Load packages
library(foreign)
library(survival)
library(rms)
library(mice)
library(glmnet)
library(lattice)
library(WGCNA)
library(Hmisc)
library(foreach)
library(doMC)
library(doParallel)
library(mstate)
library(PredictABEL)
library(ggplot2)
library(MASS)

## Load in-house functions
source("/proj/b2011036/uk.biobank/Pgm/new_function.R")

load("/proj/b2011036/uk.biobank/out8F.Rdata")

load("/proj/b2011036/uk.biobank/imputation_results/pred_selF.Rdata")

## IMPUTE ##
bdE9F <- mice2(bdE8F, m=5, predictorMatrix=pred_selF, seed=123, me="fastpmm",diagnostics=F,printFlag=T, maxit=10)

save(bdE9F,file="/proj/b2011036/uk.biobank/imputation_results/bdE9F.Rdata")

