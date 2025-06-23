#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name gwas 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas.%j.out

module load BCFtools/1.9-foss-2018a

cd /scratch/leuven/361/vsc36175/NP25_gwas_time3trials_mm80

# Remove the samples that don't have a wing measurement (not found in the phenotype.txt file)
# Shouldn't be the case with NP25
#bcftools view --threads 20 --samples-file ^remove_samples.txt gwas_imputed.vcf.gz -Oz -o gwas_imputed_clean.vcf.gz

# Checks the number of samples
bcftools query -l ../P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit.vcf.gz | wc -l

# Prepare the GWAS with PLINK
# Phenotype file should have the same FID and IID columns, and sample order should be same in vcf file 
# ! interpretation phenotype 0,1,2: for case/control status (e.g., 1=control, 2=case, 0=missing)
module load PLINK/1.9
plink --vcf ../P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit.vcf.gz --pheno ../phenotype.txt --allow-no-sex --pheno-name time3trials \
--double-id --make-bed --allow-extra-chr --out gwas_input

# To confirm that the .bed file is properly formatted
plink --bfile gwas_input --freq --allow-extra-chr --allow-no-sex

# This generates missing_check.imiss, which tells you if any individuals have missing genotype data
plink --bfile gwas_input --missing --out missing_check --allow-extra-chr --allow-no-sex

# Generates the genetic relationship matrix (GRM), also known as kinship matrix, which GEMMA requires -> gwas_input.grm.bin
plink --bfile gwas_input --make-grm-bin --out gwas_input --allow-extra-chr --allow-no-sex

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix
gemma -bfile gwas_input -gk 1 -outdir kinship_matrix -o gwas_input

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input -k kinship_matrix/gwas_input.cXX.txt -lmm 4 -outdir gemma_results -o gemma_lmm_results
