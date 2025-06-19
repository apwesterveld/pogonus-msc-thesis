#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o map_radtags.%j.out

# Submit multiple tasks as sbatch -a 1-96 script.sh

cd /scratch/leuven/361/vsc36175/

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))


# Load the programs we will use
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.18-GCC-12.3.0
module load libdeflate/1.18-GCCcore-12.3.0 # should normally load within SAMtools, but that didn't work
module load picard/2.18.23-Java-1.8.0_171
#module load minimap2/2.26-GCCcore-12.3.0

echo "================="

# Sample IDs 
#samples=(S-01/S-01 S-02/S-02 S-03/S-03 S-04/S-04 S-05/S-05 S-06/S-06 S-07/S-07 S-08/S-08 \
#S-09/S-09 S-10/S-10 S-11/S-11 S-12/S-12 S-13/S-13 S-14/S-14 S-15/S-15 S-16/S-16 \
#S-17/S-17 S-18/S-18 S-19/S-19 S-20/S-20 S-21/S-21 S-22/S-22 S-23/S-23 S-24/S-24 \
#S-25/S-25 S-26/S-26 S-27/S-27 S-28/S-28 S-29/S-29 S-30/S-30 S-31/S-31 S-32/S-32 \
#S-33/S-33 S-34/S-34 S-35/S-35 S-36/S-36 S-37/S-37 S-38/S-38 S-39/S-39 S-40/S-40 \
#S-41/S-41 S-42/S-42 S-43/S-43 S-44/S-44 S-45/S-45 S-46/S-46 S-47/S-47 S-48/S-48 \
#S-49/S-49 S-50/S-50 S-51/S-51 S-52/S-52 S-53/S-53 S-54/S-54 S-55/S-55 S-56/S-56 \
#S-57/S-57 S-58/S-58 S-59/S-59 S-60/S-60 S-61/S-61 S-62/S-62 S-63/S-63 S-64/S-64 \
#S-65/S-65 S-66/S-66 S-67/S-67 S-68/S-68 S-69/S-69 S-70/S-70 S-71/S-71 S-72/S-72 \
#S-73/S-73 S-74/S-74 S-75/S-75 S-76/S-76 S-77/S-77 S-78/S-78 S-79/S-79 S-80/S-80 \
#S-81/S-81 S-82/S-82 S-83/S-83 S-84/S-84 S-85/S-85 S-86/S-86 S-87/S-87 S-88/S-88 \
#S-89/S-89 S-90/S-90 S-91/S-91 S-92/S-92 S-93/S-93 S-94/S-94 S-95/S-95 S-96/S-96 \
#S-97/S-97 S-98/S-98 S-99/S-99 S-100/S-100 S-101/S-101 S-102/S-102 S-103/S-103 S-104/S-104 \
#S-105/S-105 S-106/S-106 S-107/S-107 S-108/S-108 S-109/S-109 S-110/S-110 S-111/S-111 S-112/S-112 \
#S-113/S-113 S-114/S-114 S-115/S-115 S-116/S-116 S-117/S-117 S-118/S-118 S-119/S-119 S-120/S-120 \
#S-121/S-121 S-122/S-122 S-123/S-123 S-124/S-124 S-125/S-125 S-126/S-126 S-127/S-127 S-128/S-128 \
#S-129/S-129 S-130/S-130 S-131/S-131 S-132/S-132 S-133/S-133 S-134/S-134 S-135/S-135 S-136/S-136 \
#S-137/S-137 S-138/S-138 S-139/S-139 S-140/S-140 S-141/S-141 S-142/S-142 S-143/S-143 S-144/S-144 \
#S-145/S-145 S-146/S-146 S-147/S-147 S-148/S-148 S-149/S-149 S-150/S-150 S-151/S-151 S-152/S-152 \
#S-153/S-153 S-154/S-154 S-155/S-155 S-156/S-156 S-157/S-157 S-158/S-158 S-159/S-159 S-160/S-160 \
#S-161/S-161 S-162/S-162 S-163/S-163 S-164/S-164 S-165/S-165 S-166/S-166 S-167/S-167 S-168/S-168 \
#S-169/S-169 S-170/S-170 S-171/S-171 S-172/S-172 S-173/S-173 S-174/S-174 S-175/S-175 S-176/S-176 \
#S-177/S-177 S-178/S-178 S-179/S-179 S-180/S-180 S-181/S-181 S-182/S-182 S-183/S-183 S-184/S-184 \
#S-185/S-185 S-186/S-186 S-187/S-187 S-188/S-188 S-189/S-189 S-190/S-190 S-191/S-191 S-192/S-192)

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
BWAout=/scratch/leuven/361/vsc36175/bams
#FILE1=/scratch/leuven/361/vsc36175/Pogonus_RAD_Nieuwpoort_Sander_18062025/$(echo "${samples[ID]}")_L1_1.fq.gz
#FILE2=/scratch/leuven/361/vsc36175/Pogonus_RAD_Nieuwpoort_Sander_18062025/$(echo "${samples[ID]}")_L1_2.fq.gz
FILE1=/scratch/leuven/361/vsc36175/Pogonus_RAD_Nieuwpoort_Sander_18062025/${samples[$ID]}/${samples[$ID]}_L1_1.fq.gz
FILE2=/scratch/leuven/361/vsc36175/Pogonus_RAD_Nieuwpoort_Sander_18062025/${samples[$ID]}/${samples[$ID]}_L1_2.fq.gz
# Alternative, map with minimap2
#minimapOUT=/scratch/leuven/361/vsc36175/bams-minimap
#file1=/scratch/leuven/361/vsc36175/$(echo "${samples[ID]}")_R1.fastq.gz
#file2=/scratch/leuven/361/vsc36175/$(echo "${samples[ID]}")_R2.fastq.gz

echo "sample ID: ${samples[$ID]}"
echo "file1 path: $FILE1"
echo "file2 path: $FILE2"
echo "reference: $REF"
echo "BWA output: $BWAout"

echo "================="

mkdir -p $BWAout

# Map reads using bwa mem
bwa mem -t 20 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/${samples[ID]}.$REFNAME.bam

# Alternative, map with minimap2
#minimap2 -ax sr -t 20 $REF $file1 $file2 | samtools view -bS - > $minimapOUT/$(echo "${samples[ID]}").$REFNAME.bam

echo "mapping finished"
echo "================="

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

echo "filtering finished"
echo "================="

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

echo "sorting finished"
echo "================="

# Remove PCR duplicates
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-INPUT $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
-OUTPUT $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-METRICS_FILE $BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt \
-ASSUME_SORTED true

echo "PCR duplicates removed"
echo "================="

# Remove intermediate files
#rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

echo "Intermediate files removed"
echo "================="

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam

echo "filtered and sorted bam indexing finished"
echo "================="
