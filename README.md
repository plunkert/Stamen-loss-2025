# Overview
This repository contains data, R scripts, bash scripts, figures, and tables associated with a QTL mapping study of short stamen number in Arabidopsis thaliana.
Cite as: Plunkert M, Issaka Salia O, Woolcock E, Perez S, Conner J. 2026. bioRxiv. Shared QTLs underlie loss of a conserved trait in replicate mapping populations of Arabidopsis thaliana.

# Scripts
1_pheno_wrangle_normalize.R takes phenotype data for RILs and parental lines (in ), wrangles it for downstream QTL mapping, and performs correlation analyses.

2_qtl_mapping_scanone_scantwo.R performs one-dimensional and two-dimensional QTL scans. Scantwo permutations computed here are used for significance thresholds in multi-QTL modeling.

3_stepwise_qtl.R performs multi-QTL modeling, producing final QTL models and associated figres and tables presented in the paper.

Directory 6_knockout_pheno_analyses contains scripts for analyzing phenotypes of TDNA insertional mutants. Phenotyping was completed in 4 separate rounds with slight differences in how raw data were recorded that require different data wrangling strategies. 

# Data
