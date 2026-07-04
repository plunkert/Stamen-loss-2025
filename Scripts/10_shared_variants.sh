#!/bin/bash --login
#SBATCH --time=5:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80GB
#SBATCH --job-name=shared_vars
#SBATCH --output=shared_vars.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=plunkert@msu.edu

# Load required packages
module purge
module load BCFtools/1.22-GCC-13.3.0 HTSlib/1.22-GCC-13.3.0

# change to directory where VCFs are stored
cd /mnt/gs21/scratch/plunkert/At_parent_reads/output/

mkdir compressed_vcfs
cp *.hc.vcf* ./compressed_vcfs
cd ./compressed_vcfs

# compress vcfs
bgzip -k Belm-12.hc.vcf
bgzip -k Kas-1.hc.vcf
bgzip -k Roda-47.hc.vcf
bgzip -k Tsu-1.hc.vcf

# index bgzipped vcfs
bcftools index Belm-12.hc.vcf.gz
bcftools index Kas-1.hc.vcf.gz
bcftools index Roda-47.hc.vcf.gz
bcftools index Tsu-1.hc.vcf.gz

module purge
module load VCFtools/0.1.16-GCC-12.3.0
vcftools --gzvcf Belm-12.hc.vcf.gz --minQ 30 --min-meanDP 10 --max-meanDP 250 --recode --stdout | gzip -c > Belm-12.filt.vcf.gz
vcftools --gzvcf Kas-1.hc.vcf.gz --minQ 30 --min-meanDP 10 --max-meanDP 250 --recode --stdout | gzip -c > Kas-1.filt.vcf.gz
vcftools --gzvcf Roda-47.hc.vcf.gz --minQ 30 --min-meanDP 10 --max-meanDP 250 --recode --stdout | gzip -c > Roda-47.filt.vcf.gz
vcftools --gzvcf Tsu-1.hc.vcf.gz --minQ 30 --min-meanDP 10 --max-meanDP 250 --recode --stdout | gzip -c > Tsu-1.filt.vcf.gz

module purge
module load BCFtools/1.22-GCC-13.3.0 HTSlib/1.22-GCC-13.3.0

# homozygous variants only
bcftools view Belm-12.filt.vcf.gz --genotype hom | bgzip -k > Belm-12.filt.hom.gz
bcftools view Kas-1.filt.vcf.gz --genotype hom | bgzip -k > Kas-1.filt.hom.gz
bcftools view Roda-47.filt.vcf.gz --genotype hom | bgzip -k > Roda-47.filt.hom.gz
bcftools view Tsu-1.filt.vcf.gz --genotype hom | bgzip -k > Tsu-1.filt.hom.gz

bcftools index Belm-12.filt.hom.gz
bcftools index Kas-1.filt.hom.gz
bcftools index Roda-47.filt.hom.gz
bcftools index Tsu-1.filt.hom.gz


# get variants shared by Belm and Kas but not Roda and Tsu. Have to be the same variant, not just in the same place
bcftools isec -n~1100 -c none Belm-12.filt.hom.vcf.gz Kas-1.filt.hom.vcf.gz Roda-47.filt.hom.vcf.gz Tsu-1.filt.hom.vcf.gz > stamen_shared_vars.vcf
