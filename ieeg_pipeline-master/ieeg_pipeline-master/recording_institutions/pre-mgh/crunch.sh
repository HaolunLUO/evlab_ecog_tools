#!/bin/bash
#SBATCH --job-name=ieeg_preprocessing
#SBATCH --time=03:00:00
#SBATCH --mem=128GB
#SBATCH --partition=ou_bcs_normal
#SBATCH --account=mit_general
#SBATCH --output="logs/ieeg_preprocessing_%j.out"
#SBATCH -c 4

echo $1 # subject file
echo $2 # order file
echo $3 # experiment name
echo $SLURMD_NODENAME

subject=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${1} )
order=$( awk "NR==$SLURM_ARRAY_TASK_ID" ${2} )

mkdir -p logs
module load matlab/matlab-2025b
matlab -nodisplay -r "crunch('$subject','${3}','fromScratch',true,'doneVisualInspection',true,'order','$order','isPlotVisible',false);"
