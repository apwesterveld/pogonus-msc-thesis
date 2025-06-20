#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name filter 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o filter_vcf.%j.out

# 20250620 VCFtools (0.1.16)
source /user/leuven/361/vsc36175/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
module load tabix/0.2.6-GCCcore-6.4.0

cd /scratch/leuven/361/vsc36175/

# Check missingness per individual
vcftools --gzvcf P_chalceus_NP25_BarSW_merged.vcf.gz --missing-indv

# Filtering
vcftools --gzvcf P_chalceus_NP25_BarSW_merged.vcf.gz \
--max-missing 1 --minQ 30 --maf 0.05 --remove-indels --recode --stdout | bgzip > P_chalceus_NP25_BarSW_merged_filtered.vcf.gz
# max-missing 0.9: At least 90% of individuals must have a called genotype. First try with 100%
# minQ 30: minimum quality score of 30
# maf 0.05: minor allele frequency (MAF) threshold of 0.05
# remove-indels: exclude indels, keeping only SNPs

# Split the multiallelic SNPs into multiple biallelic ones
module load BCFtools/1.9-foss-2018a
bcftools norm -m -any -o P_chalceus_NP25_BarSW_merged_filtered_multiSplit.vcf.gz -Oz P_chalceus_NP25_BarSW_merged_filtered.vcf.gz

# Index VCF file
tabix -p vcf P_chalceus_NP25_BarSW_merged_filtered_multiSplit.vcf.gz
