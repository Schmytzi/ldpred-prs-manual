#!Rscript

# Script does not accept command line paramters. Please customize it before running
# I recommend pasting it step-by-step into an interactive R session.
# You might want to use some intermediate results
# Use AT LEAST 8 cores on Bianca. The combined dataset is huge

library(data.table)
library(magrittr)
library(stringr)

# Set to TRUE for quantitative phenotypes
is_binary <- FALSE

get_gwas <- function(file_pattern){
    list.files(path="gwas", pattern=file_pattern, full.names=TRUE) %>%
        lapply(fread) %>%
        rbindlist()
}

convert_to_ldpred <- function(data, effect_column) {
    data[, .(
        chr = CHR,
        pos = POS.x,
        ref = A1.x,
        # A2 is a comma-separated list of all non-A1 alleles
        # If the input bim hasn't been deduplicated, this can lead to invalid output 
        alt = A2,
        # reffrq is not the reference's freq, but the reference minor allele frequency
        reffrq = MAF,
        info = INFO,
        rs = rsIds,
        pval = P,
        effalt = effect_column
    )]
}

# Load pre-computed SNP metadata
mfi <- list.files(path="mfi", pattern="*mfi_chr*", full.names=TRUE) %>% 
    lapply(fread, header=FALSE) %>% 
    rbindlist()
names(mfi) <- c("SNP", "ID", "POS", "A1", "A2", "MAF", "MA", "INFO")

# Load GWAS summary
gwas <- get_gwas(if (is_binary) "*logistic" else "*linear")
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

ldpred <- convert_to_ldpred(merged, if (is_binary) merged$OR else merged$BETA)
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