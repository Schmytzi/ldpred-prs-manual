#!/bin/bash -l

#SBATCH -A sens2017538
#SBATCH -J coord
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 7:00:00

module load bioinfo-tools LDpred/1.0.8
ldpred coord \
  --gf reference_genotypes/merged \
  --ssf gwas/LDpred_format.txt \
  --ssf-format CUSTOM \
  --chr chr \
  --pos pos \
  --info info \
  --rs rs \
  --A1 a1 \
  --A2 a2 \
  --reffreq reffreq \
  --pval pval \
  --eff eff \
  --ncol nsamples \
  --beta \
  --out coordination/coordination.coord
