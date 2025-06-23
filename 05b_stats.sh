#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name stats 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o stats.%j.out

module load BCFtools/1.9-foss-2018a

source /user/leuven/361/vsc36175/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

cd /scratch/leuven/361/vsc36175/

# Gives variant count, transition/transversion ratio, missing genotype rates, etc..
bcftools stats P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz

# Gives a brief summary of stats
bcftools stats P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz | grep -E "SN|TSTV"

# Counts the number of variants (excluding the header).
bcftools view -H P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz | wc -l

# Checks missingness
vcftools --gzvcf P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz --missing-indv

# Checks allele frequency distribution 
vcftools --gzvcf P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz --freq

# Checks allele density
vcftools --gzvcf P_chalceus_NP25_BarSW_merged_filtered_mm80_thinned5k_multiSplit.vcf.gz --SNPdensity 1000000 --out snp_density
