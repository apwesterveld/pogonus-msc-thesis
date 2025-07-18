#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name snps 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out

# Submit multiple tasks as: sbatch -a 1-96 script.sh

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index  starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BCFtools/1.9-foss-2018a

#for local BCFtools:
#module load Python/3.7.0-foss-2018a
#export BCFTOOLS_PLUGINS=/data/leuven/361/vsc36175/bcftools/plugins

echo "================="

chrom=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10)
names=(01 02 03 04 05 06 07 08 09 10)

# Sample IDs
samples=(S-01 S-02 S-03 S-04 S-05 S-06 S-07 S-08 \
S-09 S-10 S-11 S-12 S-13 S-14 S-15 S-16 \
S-17 S-18 S-19 S-20 S-21 S-22 S-23 S-24 \
S-25 S-26 S-27 S-28 S-29 S-30 S-31 S-32 \
S-33 S-34 S-35 S-36 S-37 S-38 S-39 S-40 \
S-41 S-42 S-43 S-44 S-45 S-46 S-47 S-48 \
S-49 S-50 S-51 S-52 S-53 S-54 S-55 S-56 \
S-57 S-58 S-59 S-60 S-61 S-62 S-63 S-64 \
S-65 S-66 S-67 S-68 S-69 S-70 S-71 S-72 \
S-73 S-74 S-75 S-76 S-77 S-78 S-79 S-80 \
S-81 S-82 S-83 S-84 S-85 S-86 S-87 S-88 \
S-89 S-90 S-91 S-92 S-93 S-94 S-95 S-96 \
S-97 S-98 S-99 S-100 S-101 S-102 S-103 S-104 \
S-105 S-106 S-107 S-108 S-109 S-110 S-111 S-112 \
S-113 S-114 S-115 S-116 S-117 S-118 S-119 S-120 \
S-121 S-122 S-123 S-124 S-125 S-126 S-127 S-128 \
S-129 S-130 S-131 S-132 S-133 S-134 S-135 S-136 \
S-137 S-138 S-139 S-140 S-141 S-142 S-143 S-144 \
S-145 S-146 S-147 S-148 S-149 S-150 S-151 S-152 \
S-153 S-154 S-155 S-156 S-157 S-158 S-159 S-160 \
S-161 S-162 S-163 S-164 S-165 S-166 S-167 S-168 \
S-169 S-170 S-171 S-172 S-173 S-174 S-175 S-176 \
S-177 S-178 S-179 S-180 S-181 S-182 S-183 S-184 \
S-185 S-186 S-187 S-188 S-189 S-190 S-191 S-192)

# Some folder and file paths to use later
REF=/scratch/leuven/361/vsc36175/reference_genome/P_chalceus.fasta
REFNAME=BarSW

cd /scratch/leuven/361/vsc36175/bams

# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in "${samples[@]}"; do
    ALL_LIST+=" ${FILE}.${REFNAME}.filtered.sorted.bam"
done
 
command="$ALL_LIST"

echo "Reference: $REF"
echo "Chromosome: ${chrom[ID]}"
echo "Command: $command"

# Ensure chromosome ID is not empty
if [[ -z "${chrom[ID]}" ]]; then
    echo "Error: chrom[ID] is empty!"
    exit 1
fi

# Ensure BAM files exist
if [[ -z "$command" ]]; then
    echo "Error: No BAM files specified!"
    exit 1
fi

# run mpileup
bcftools mpileup -Oz --threads 20 -f $REF $command -r ${chrom[ID]} | \
bcftools call -m -Oz -o /scratch/leuven/361/vsc36175/P_chalceus_NP25_$REFNAME.chr_${names[ID]}.vcf.gz
