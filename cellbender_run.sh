#!/bin/bash 
#BSUB -J run_4
#BSUB -P acc_psychgen 
#BSUB -q gpuexpress
#BSUB -n 1 
#BSUB -W 12:00
#BSUB -R rusage[mem=20G]
#BSUB -R affinity[core(32)]
#BSUB -gpu num=1
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml purge
unset PYTHONPATH
ml anaconda3/2021.5
source activate /sc/arion/work/busson02/LBP/envs/CellBender
unset PYTHONPATH


mkdir PT-292-blood-L

cellbender remove-background  --input /sc/arion/projects/psychgen/lbp/data/RAW/rna/sc/chromium/TD005881/TD005881_AlexCharney/R01AG069976AIM2SC1761/outs/raw_feature_bc_matrix.h5 --output PT-290-blood-L/PT-290-L-CellBender.h5  --cuda