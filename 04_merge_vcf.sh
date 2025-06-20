#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name merge 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_vcfs.%j.out

module load BCFtools/1.9-foss-2018a

# For local bcftools install:
#module load Python/3.7.0-foss-2018a
#export BCFTOOLS_PLUGINS=/data/leuven/361/vsc36175/bcftools/plugins

cd /scratch/leuven/361/vsc36175/

bcftools concat -Oz -o P_chalceus_NP25_BarSW_merged.vcf.gz *.vcf.gz
