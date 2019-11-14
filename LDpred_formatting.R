#!Rscript

# Script does not accept command line paramters. Please customize it before running
# I recommend pasting it step-by-step into an interactive R session.
# You might want to use some intermediate results
# Use AT LEAST 8 cores on Bianca. The combined dataset is huge

library(data.table)
library(magrittr)
library(stringr)

# Load pre-computed SNP metadata
mfi <- list.files(path="mfi", pattern="*mfi_chr*", full.names=TRUE) %>% 
    lapply(fread, header=FALSE) %>% 
    rbindlist()
names(mfi) <- c("SNP", "ID", "POS", "A1", "A2", "MAF", "MA", "INFO")

# Load GWAS summary
gwas <- list.files(path="gwas", pattern="*linear", full.names=TRUE) %>%
    lapply(fread) %>%
    rbindlist()
# Discard all rows without results
gwas <- gwas[complete.cases(gwas)]
# Sometimes, we must force R to treat the P values as number
gwas[, P := as.numeric(P)]
# #CHROM is difficult to work with as a name
setnames(gwas, "#CHROM", "CHR")

# Use this if you haven't used deduplicated bim files
  merged <- merge(mfi, gwas, by.x="ID", by.y="ID")
  # Column ID contains rs IDs
  rsIds <- merged$ID
# Use this if you have used deduplicated bim files
  matches <- str_match(gwas$ID, "(rs\\d+)_\\d")
  gwas[, oldID := ifelse(is.na(matches[,2]), ID, matches[,2])]
  merged <- merge(mfi, gwas, by.x="ID", by.y="oldID")
  # Merged on old IDs, so the deduplicated rsIDs are in column ID.y
  rsIds <- merged$ID.y

ldpred <- merged[, .(
    chr = CHR,
    pos = POS.x,
    ref = REF,
    alt = ALT,
    # If the reference is the minor allele -> MAF = reffrq, else MAF refers to the alt's freq
    reffrq = ifelse(A1.y == REF, MAF, 1 - MAF),
    info = INFO,
    rs = rsIds,
    pval = P,
    effalt = BETA
)]
# Ensure proper ordering
setkey(ldpred, chr, pos)
# LDpred expects all chromosome specifications to start with "chr"
ldpred[, chr := paste0("chr", chr)]

# FOR LDPRED < 1.0.0 ONLY:
# LDpred version 0.9.9 cannot handle P values that have been rounded to 0.
# Therefore, replace all 0s with values close to the minimum possible.
# Uncomment the following line if you're using an old version of LDpred
# ldpred[pval == 0, pval := 2.23e-308]

write.table(ldpred, "gwas/LDpred_format.txt", quote=FALSE, row.names=FALSE)