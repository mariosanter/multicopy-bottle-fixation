#!/bin/bash
#SBATCH --job-name=selection
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=01:00:00
#SBATCH --output=batch-%A_%a.out
#SBATCH --error=batch-%A_%a.err
#SBATCH --partition=base
#SBATCH --array=0-2400
OMP_NUM_THREADS=1
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate bottlenecks2023

python - << EOF
import module
import json
p=json.load(open( "selection-SI_input/inputfile_${SLURM_ARRAY_TASK_ID}", 'r'))
del p["repindex"]
sim=module.stochbottleSim(**p)
EOF
