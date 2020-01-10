#!/bin/bash -l

#SBATCH -A sens2017538
#SBATCH -J coord
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 7:00:00
#SBATCH -C fat

# DO NOT FORGET to provide a valid value for N
module load bioinfo-tools LDpred/1.0.8
ldpred coord \
  --gf reference_genotypes/merged \
  --ssf gwas/LDpred_format.txt \
  --N 0 \
  --ssf-format CUSTOM \
  --chr chr \
  --pos pos \
  --a1 a1 \
  --a2 a2 \
  --reffreq reffrq \
  --info info \
  --rs rs \
  --pval pval \
  --eff effalt \
  --out coordination/coordination.coord
