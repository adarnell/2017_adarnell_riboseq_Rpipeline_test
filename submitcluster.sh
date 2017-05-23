#!/bin/bash -l

#SBATCH --mem=8000
#SBATCH -t 08:00:00  
#SBATCH -J alicia

make -f Makefile.align all SAMPLE=$1
