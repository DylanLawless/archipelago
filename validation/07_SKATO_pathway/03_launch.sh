#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 200G
#SBATCH --time 32:00:00
#SBATCH --job-name=archi
#SBATCH --output=./log_launch_%a_%A_%J.out
#SBATCH --error=./log_launch_%a_%A_%J.err
#SBATCH --partition=dynamic-64cores-256g-1gpu-A100-40g

set -e
module load r/4.3.1
echo "START AT $(date)"

Rscript 04_run_skato.R
echo "FINISH AT $(date)"
