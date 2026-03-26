#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1

#
chromosome=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/chromosomes_file`
inputdir=/path/to/CG_reports

for file in ${inputdir}/*.CG_report.txt; do
    grep -w "${chromosome}" "$file" > "${inputdir}/${chromosome}/$(basename "${file}" .CG_report.txt)_${chromosome}_CG_report.txt"
done
