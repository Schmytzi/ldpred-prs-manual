#!/bin/bash -l
#SBATCH -A sens2017538
#SBATCH -J Scoring
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 8:00:00
#SBATCH -a 0-175
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.schmitz@igp.uu.se

module load bioinfo-tools plink/1.90b4.9 
cd /home/danies/shbg/plink2/males/ld_estimation/weights

let "chr = $SLURM_ARRAY_TASK_ID % 22 + 1"
let "col = $SLURM_ARRAY_TASK_ID / 22 + 6"
col_name=$(head -n 1 all_weights.txt | cut -d ' ' -f $col)

echo "Scoring chromosome $chr and column $col_name"

mkdir -p out/col$col

plink \
  --bed ~/shbg/plink2/dosages/chr${chr}.bed \
  --bim ~/shbg/plink2/dosages/chr${chr}.dedup.bim \
  --fam /proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3/ukb41143.fam \
  --score all_weights.txt 3 5 $col header sum \
  --extract ../info_snps.txt \
  --filter-males \
  --keep /proj/sens2017538/nobackup/UKBB_IMP_DOSAGE_V3/IDs_to_indlude_NonRelCauc.ukb41143.txt \
  --out out/col$col/scores_chr${chr}
