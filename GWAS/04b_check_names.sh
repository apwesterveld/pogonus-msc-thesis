#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name check_names 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o check_names.%j.out

module load BCFtools/1.9-foss-2018a

cd /scratch/leuven/361/vsc36175/

bcftools query -l P_chalceus_NP25_BarSW_merged.vcf.gz > vcf_samples.txt
cut -f2 phenotype_modified.txt  > phenotype_samples.txt
comm -3 <(sort vcf_samples.txt) <(sort phenotype_samples.txt)
