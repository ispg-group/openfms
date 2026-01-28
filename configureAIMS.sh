#!/bin/bash -l
#SBATCH --job-name=ConfigAIMS
#SBATCH --partition=rooster,chem          
#SBATCH --nodes=1
#SBATCH --ntasks=1                  
#SBATCH --time=04:00:00
#SBATCH --output=ConfigAIMS.%j.o
#SBATCH --error=ConfigAIMS.%j.e
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=REQUEUE
##SBATCH --nodelist=gpu193
#SBATCH --mem=200GB


./configure.sh
make -j
make unittest
make test 

