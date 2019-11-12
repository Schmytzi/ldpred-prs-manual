#!Rscript

library(data.table)

bed <- list.files(path="reference_genotypes", pattern="bed")
bim <- list.files(path="reference_genotypes", pattern="bim")
# plink generates new fam files for each chromosome
fam <- list.files(path="reference_genotypes", pattern="fam")

files <- data.table(bed, bim, fam)
# Omit line with chr1 b.c. that will be the main data set for plink
files <- files[!grepl("chr1.bed", bed)]
write.table(
    files,
    "reference_genotypes/merge_list.txt",
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE)