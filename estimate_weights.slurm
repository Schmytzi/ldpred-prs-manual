#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -J ldpred
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00

module load bioinfo-tools LDpred/1.0.8
ldpred gibbs \
  --ldr 0 \
  --ldf ld/ld_pattern \
  --N 0 \
  --coord coordination/coordination.coord \
  --out weights/weights
