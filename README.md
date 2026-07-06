# Overview
This repository contains data, R scripts, bash scripts, figures, and tables associated with a QTL mapping study of short stamen number in Arabidopsis thaliana.
Cite as: Plunkert M, Issaka Salia O, Woolcock E, Perez S, Conner J. 2026. bioRxiv. Shared QTLs underlie loss of a conserved trait in replicate mapping populations of Arabidopsis thaliana.

# Scripts
1_pheno_wrangle_normalize.R takes phenotype data for RILs and parental lines (in ), wrangles it for downstream QTL mapping, and performs correlation analyses.

2_percent_parental_contribution.R estimates the percent of the genome of each RIL that came from each parent, and calculates partial correlations between phenotypes that are corrected for parental genome contribution.

3_qtl_mapping_scanone_scantwo.R performs one-dimensional and two-dimensional QTL scans. Scantwo permutations computed here are used for significance thresholds in multi-QTL modeling.

4_stepwise_qtl.R performs multi-QTL modeling, producing final QTL models and associated figres and tables presented in the paper.

5_investigate_linked_QTL.R examines pairs of putative linked QTLs on chromosome 1 and 5. Linked QTLs can be dubious, and so this analysis is not emphasized in the paper.

Directory 6_knockout_pheno_analyses contains scripts for analyzing phenotypes of TDNA insertional mutants. Phenotyping was completed in 4 separate rounds with slight differences in how raw data were recorded that require different data wrangling strategies, hence the separate scripts here. 

7_genomic_region_figure.R produces a figure showing the genomic region from which candidate genes were selected.

8_map_reads_call_vars.sh takes FASTQ data, maps to the TAIR10 reference genome, and calls variants. Input data are on SRA (BioProject PRJNA1468757).

9_shared_variants.sh performs some VCF filtering and finds homozygous variants that are found in stamen loss but not stamen retaining parents.

10_shared_vars_snpEff.sh uses snpEff to predict effects of the shared variants on gene function.

# Data

Directory belm_roda_pheno_geno contains phenotype and genotype data for Belm-12 x Roda-47 RILs.
- Chamber_FloweringTime.xlsx contains flowering time data from Dittmar et al., 2014, Flowering time QTL in natural populations of Arabidopsis thaliana and implications for their adaptive value, in Molecular Ecology.
- IT SW RIL genotypes.csv contains RIL genotypes.
- RIL mean ShortStamenNo.csv contains mean short stamen number for RILs.
- Royer final QTL data for analysis.csv contains data from Royer et al., 2016, Incomplete loss of a conserved trait: function, latitudinal cline, and genetic constraints, Evolution.
- SW IT OvuleNo.xlsx contains ovule number data for Belm-12 x Roda-47 RILs and parental lines
- SW_IT_RILparents_forVariance.csv contains stamen number data for Belm-12 and Roda-47 parental lines.
- wrangled_belm_roda_phenos.csv contains RIL phenotypes following transformation and other data wrangling. Input file for QTL mapping.

Directory tsu_kas_pheno_geno contains phenotype and genotype data for Tsu-1 x Kas-1 RILs.
- Tsu_Kas_Rqtl_format.csv contains genotype and phenotype data in Rqtl format.
- TsuKas_RIL_ovule_counts.xlsx contains ovule count data for Tsu-1 x Kas-1 RIL
- Tsu-Kas RILs.xls contains raw RIL phenotype and parent data.
- TKrilsGenos_map55_pheno.csv
- TK_RILs_IDs.csv contains the multiple IDs used for Tsu-1 x Kas-1 RILs on ABRC and in various genetic maps.
- TsuKas_parent_ovule_counts_MLPentry.xlsx

Directory knockout_phenos contains phenotype data for all 4 rounds of TDNA insertional mutant phenotyping, numbered 1-4.






