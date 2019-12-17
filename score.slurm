#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -J Scoring
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 8:00:00
#SBATCH -a 0-175

module load bioinfo-tools plink/1.90b4.9 

let "chr = $SLURM_ARRAY_TASK_ID % 22 + 1"
let "col = $SLURM_ARRAY_TASK_ID / 22 + 6"
col_name=$(head -n 1 all_weights.txt | cut -d ' ' -f $col)

echo "Scoring chromosome $chr and column $col_name"

mkdir -p scores/col$col

plink \
  --bed genotypes/chr${chr}.bed \
  --bim genotypes/chr${chr}.dedup.bim \
  --fam sample_inf/sample.fam \
  --score weights/all_weights.txt 3 5 $col header sum \
  --keep sample_info/validation.ids \
  --out scores/col$col/scores_chr${chr}
