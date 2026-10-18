#!/bin/bash
#BSUB -J train_L2R          # job name
#BSUB -W 60                 # wall clock time limit (miniutes)
#BSUB -n 3                  # number of CPU cores required
#BSUB -R "rusage[mem=5GB]"  # memory

#BSUB -q ccee               # name of queue for submission
#BSUB -o result_%J.out      # out put file
#BSUB -e result_%J.err      # error file

source ~/.bashrc

echo "activating conda environment........"
conda activate $ENV_PATH

# run script
python train_L2R.py

conda deactivate

echo "Job finished" 