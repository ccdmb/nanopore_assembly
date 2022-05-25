#!/bin/bash --login

#SBATCH --partition=workq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --export=NONE

module load nextflow
module load java

srun --export=all nextflow run -profile singularity,zeus -resume ./main.nf --nanoporeReads '$*'
