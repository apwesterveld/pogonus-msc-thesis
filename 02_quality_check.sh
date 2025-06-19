#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name qc
#SBATCH --nodes=1 
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o bams_qc.%j.out

cd /scratch/leuven/361/vsc36175/bams

module load SAMtools/1.9-GCC-6.4.0-2.28

for bam in *.filtered.sorted.np.bam; do 
    echo "$bam"; 
    samtools flagstat "$bam"; 
done > mapping_stats.txt

# Output file for errors
ERROR_LOG="bam_quickcheck_errors.txt"

# Empty the error log file before running
> "$ERROR_LOG"

# Loop through all filtered.sorted.bam files
for bam in *.filtered.sorted.bam; do
    if ! samtools quickcheck "$bam"; then
        echo "$bam" >> "$ERROR_LOG"
    fi
done

echo "Check complete. Files with errors are listed in: $ERROR_LOG"

for bam in *.filtered.sorted.nd.bam; do 
    echo -n "$(basename "$bam" .bam): "; 
    samtools depth -a "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print "No coverage"}'; 
done > average_depths.txt
