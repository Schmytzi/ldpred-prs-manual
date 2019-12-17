#!/bin/bash -l

#SBATCH -A sens2017538
#SBATCH -J LD-Extract-Individuals
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 3:00:00
#SBATCH -a 1-22

chr=$SLURM_ARRAY_TASK_ID

module load bioinfo-tools plink2/2.00-alpha-2-20190429

# If you don't want to filter by INFO score, remove the --extract parameter
plink2 \
  --pgen genotypes/chr${chr}.pgen \
  --psam sample_info/samples.fam \
  --pvar genotypes/chr${chr}.dedup.bim \
  --keep sample_info/estimation.ids \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-20 \
  --extract info_snps.txt \
  --make-bed \
  --out reference_genotypes/chr${chr}