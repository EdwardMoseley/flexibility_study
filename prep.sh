#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --error=submit.err
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4 #48
#SBATCH --partition=compsci

set FILENAME = $1
set ARG = "$FILENAME"
set CMD = "python3 prep.py $ARG"
eval $CMD
echo "Completed!"
