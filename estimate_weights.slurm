#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -J ldpred
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.schmitz@igp.uu.se

module load bioinfo-tools LDpred/1.0.8
ldpred gibbs \
  --ld_radius 0 \
  --local_ld_file_prefix UKB_LD \
  --N 0 \
  --coord coordination/coordination.coord \
  --out weights/weights
