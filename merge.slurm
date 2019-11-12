#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -J Merge
#SBATCH -p node
#SBATCH -N 1
#SBATCH -C fat
#SBATCH -t 17:00:00

module load bioinfo-tools plink/1.90b3n

# We require a fat node and tell plink how much memory to use
# Otherwise we'd run out of memory
plink \
  --bfile reference_genotypes/chr1 \
  --merge-list reference_genotypes/merge_list.txt \
  --make-bed \
  --memory 220000 \
  --out reference_genotypes/merged
