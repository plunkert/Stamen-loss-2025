#!/bin/bash

cd ~/Documents/GitHub/Stamen-loss-2025/

# path to snpEff program files
snpeff=~/Simple-1.8.1-master/programs/snpEff/

java -jar $snpeff/snpEff.jar Arabidopsis_thaliana -s ./Results/snpEff_summary_test.html ./Results/stamen_shared_vars.vcf > ./Results/stamen_shared_vars_ann.vcf