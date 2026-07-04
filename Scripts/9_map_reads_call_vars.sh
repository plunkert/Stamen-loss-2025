#!/bin/bash --login
#SBATCH --time=19:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80GB
#SBATCH --job-name=map_At_parent_reads
#SBATCH --output=map_At_parent.SLURMout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=plunkert@msu.edu

# Load required packages
module purge
#module load Miniforge3/25.11.0-1
#conda activate mapIlluminaReads2
#module load SAMtools/1.17-GCC-12.2.0
#module load Java/17.0.6

#set variables
# change to path where TAIR10 reference genome is stored
ref=/mnt/gs21/scratch/plunkert/At_parent_reads/refs/TAIR10_chr_all.fasta
threads=10

# Trim reads using fastp
fastp -i parent_raw_reads/Belm-12_gDNA_S3_L003_R1_001.fastq.gz -I parent_raw_reads/Belm-12_gDNA_S3_L003_R2_001.fastq.gz -o output/Belm-12_trimmed_1.fastq.gz -O output/Belm-12_trimmed_2.fastq.gz
fastp -i parent_raw_reads/Kas-1_gDNA_S5_L003_R1_001.fastq.gz -I parent_raw_reads/Kas-1_gDNA_S5_L003_R2_001.fastq.gz -o output/Kas-1_trimmed_1_fastq.gz -O output/Kas-1_trimmed_2.fastq.gz
fastp -i parent_raw_reads/Roda-47_gDNA_S4_L003_R1_001.fastq.gz -I parent_raw_reads/Roda-47_gDNA_S4_L003_R2_001.fastq.gz -o output/Roda-47_trimmed_1_fastq.gz -O output/Roda-47_trimmed_2.fastq.gz
fastp -i parent_raw_reads/Tsu-1_gDNA_S6_L003_R1_001.fastq.gz -I parent_raw_reads/Tsu-1_gDNA_S6_L003_R2_001.fastq.gz -o output/Tsu-1_trimmed_1.fastq.gz -O output/Tsu-1_trimmed_2.fastq.gz


####creating reference files#### uncomment if ref genome not already indexed
# creating .fai file
#samtools faidx refs/TAIR10_chr_all.fas
# creating bwa-mem2 index files
#bwa-mem2 index refs/TAIR10_chr_all.fas
#mv TAIR10_chr* refs/

#generating dict file for GATK # comment out after first use
#java -Xmx2g -jar $PICARD/picard.jar CreateSequenceDictionary R=refs/TAIR10_chr_all.fas O=refs/TAIR10_chr_all.dict

# map reads to reference
bwa-mem2 mem -t ${threads} ${ref} output/Belm-12_trimmed* > output/Belm-12.sam
bwa-mem2 mem -t ${threads} ${ref} output/Kas-1_trimmed* > output/Kas-1.sam
bwa-mem2 mem -t ${threads} ${ref} output/Roda-47_trimmed* > output/Roda-47.sam
bwa-mem2 mem -t ${threads} ${ref} output/Tsu-1_trimmed* > output/Tsu-1.sam

wait

samtools view -bSh output/Belm-12.sam > output/Belm-12.bam
samtools view -bSh output/Kas-1.sam > output/Kas-1.bam
samtools view -bSh output/Roda-47.sam > output/Roda-47.bam
samtools view -bSh output/Tsu-1.sam > output/Tsu-1.bam

samtools fixmate output/Belm-12.bam output/Belm-12.fix.bam
samtools fixmate output/Kas-1.bam output/Kas-1.fix.bam
samtools fixmate output/Roda-47.bam output/Roda-47.fix.bam
samtools fixmate output/Tsu-1.bam output/Tsu-1.fix.bam

wait

samtools sort -o output/Belm-12.sort.bam output/Belm-12.fix.bam
samtools sort -o output/Kas-1.sort.bam output/Kas-1.fix.bam
samtools sort -o output/Roda-47.sort.bam output/Roda-47.fix.bam
samtools sort -o output/Tsu-1.sort.bam output/Tsu-1.fix.bam

wait


#this part is just to add header for further gatk tools. Trying it on .sort.bam files before MarkDuplicates step because having issues with MarkDuplicates
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Belm-12.sort.bam O=output/Belm-12.sort.rg.bam RGLB=Belm-12 RGPL=illumina RGSM=Belm-12 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Kas-1.sort.bam O=output/Kas-1.sort.rg.bam RGLB=Kas-1 RGPL=illumina RGSM=Kas-1 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Roda-47.sort.bam O=output/Roda-47.sort.rg.bam RGLB=Roda-47 RGPL=illumina RGSM=Roda-47 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Tsu-1.sort.bam O=output/Tsu-1.sort.rg.bam RGLB=Tsu-1 RGPL=illumina RGSM=Tsu-1 RGPU=run1 SORT_ORDER=coordinate
wait

java -Xmx2g -jar $PICARD/picard.jar MarkDuplicates I=output/Belm-12.sort.rg.bam O=output/Belm-12.sort.md.rg.bam METRICS_FILE=output/Belm-12.matrics.txt ASSUME_SORTED=true
java -Xmx2g -jar $PICARD/picard.jar MarkDuplicates I=output/Kas-1.sort.rg.bam O=output/Kas-1.sort.md.rg.bam METRICS_FILE=output/Kas-1.matrics.txt ASSUME_SORTED=true
java -Xmx2g -jar $PICARD/picard.jar MarkDuplicates I=output/Roda-47.sort.rg.bam O=output/Roda-47.sort.md.rg.bam METRICS_FILE=output/Roda-47.matrics.txt ASSUME_SORTED=true
java -Xmx2g -jar $PICARD/picard.jar MarkDuplicates I=output/Tsu-1.sort.rg.bam O=output/Tsu-1.sort.md.rg.bam METRICS_FILE=output/Tsu-1.matrics.txt ASSUME_SORTED=true

wait

#this part is just to add header for further gatk tools
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Belm-12.sort.md.bam O=output/Belm-12.sort.md.rg.bam RGLB=Belm-12 RGPL=illumina RGSM=Belm-12 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Kas-1.sort.md.bam O=output/Kas-1.sort.md.rg.bam RGLB=Kas-1 RGPL=illumina RGSM=Kas-1 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Roda-47.sort.md.bam O=output/Roda-47.sort.md.rg.bam RGLB=Roda-47 RGPL=illumina RGSM=Roda-47 RGPU=run1 SORT_ORDER=coordinate
java -Xmx2g -jar $PICARD/picard.jar AddOrReplaceReadGroups I=output/Tsu-1.sort.md.bam O=output/Tsu-1.sort.md.rg.bam RGLB=Tsu-1 RGPL=illumina RGSM=Tsu-1 RGPU=run1 SORT_ORDER=coordinate

wait


# Filter read pairs where one mate mapped and not the other
samtools view -F 284 output/Belm-12.sort.md.rg.bam -o output/Belm-12.sort.md.rg.filt.bam
samtools view -F 284 output/Kas-1.sort.md.rg.bam -o output/Kas-1.sort.md.rg.filt.bam
samtools view -F 284 output/Roda-47.sort.md.rg.bam -o output/Roda-47.sort.md.rg.filt.bam
samtools view -F 284 output/Tsu-1.sort.md.rg.bam -o output/Tsu-1.sort.md.rg.filt.bam

java -Xmx2g -jar $PICARD/picard.jar BuildBamIndex INPUT=output/Belm-12.sort.md.rg.filt.bam
java -Xmx2g -jar $PICARD/picard.jar BuildBamIndex INPUT=output/Kas-1.sort.md.rg.filt.bam
java -Xmx2g -jar $PICARD/picard.jar BuildBamIndex INPUT=output/Roda-47.sort.md.rg.filt.bam
java -Xmx2g -jar $PICARD/picard.jar BuildBamIndex INPUT=output/Tsu-1.sort.md.rg.filt.bam


# Index bam before variant calling with GATK
samtools index output/Belm-12.sort.md.rg.filt.bam
samtools index output/Kas-1.sort.md.rg.filt.bam
samtools index output/Roda-47.sort.md.rg.filt.bam
samtools index output/Tsu-1.sort.md.rg.filt.bam

# Deactivate conda environment
conda deactivate

# Load GATK for variant calling
module purge
module load GATK/4.5.0.0-GCCcore-12.3.0-Java-17

# Call variants

gatk --java-options "-Xmx40G" HaplotypeCaller --native-pair-hmm-threads 10 -R ${ref} -I output/Belm-12.sort.md.rg.filt.bam -I output/Kas-1.sort.md.rg.filt.bam -I output/Roda-47.sort.md.rg.filt.bam -I output/Tsu-1.sort.md.rg.filt.bam --min-base-quality-score 20 -O output/RIL_parents.hc.vcf

gatk --java-options "-Xmx40G" HaplotypeCaller --native-pair-hmm-threads 10 -R ${ref} -I output/Belm-12.sort.md.rg.filt.bam --min-base-quality-score 20 -O output/Belm-12.hc.vcf
gatk --java-options "-Xmx40G" HaplotypeCaller --native-pair-hmm-threads 10 -R ${ref} -I output/Kas-1.sort.md.rg.filt.bam --min-base-quality-score 20 -O output/Kas-1.hc.vcf
gatk --java-options "-Xmx40G" HaplotypeCaller --native-pair-hmm-threads 10 -R ${ref} -I output/Roda-47.sort.md.rg.filt.bam --min-base-quality-score 20 -O output/Roda-47.hc.vcf
gatk --java-options "-Xmx40G" HaplotypeCaller --native-pair-hmm-threads 10 -R ${ref} -I output/Tsu-1.sort.md.rg.filt.bam --min-base-quality-score 20 -O output/Tsu-1.hc.vcf

# Print resource information
scontrol show job $SLURM_JOB_ID
js -j $SLURM_JOB_ID

