#!/bin/bash
#SBATCH --qos=batch
#SBATCH --partition=compute
#SBATCH --job-name=scRATE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=23:59:59
#SBATCH --array=1-6

module load singularity

# Run your R script that uses scRATE singularity container
ARRAY_ID=`printf %05d $SLURM_ARRAY_TASK_ID`
singularity run --app Rscript ${CONTAINER} ${RFILE} _chunk.${ARRAY_ID} _scrate.${ARRAY_ID} ${CORES} ${SEED}
