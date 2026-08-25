#!/usr/bin/env Rscript
# =====================================================================
# p5_summarize.R - Summarize all population private SNP statistics (read private_out)
# Usage: Rscript p5_summarize.R <outdir> [popfile]
# Output: summary.tsv / private_multipop.tsv / summary.png
# =====================================================================
args    <- commandArgs(trailingOnly = TRUE)
OUTDIR  <- if (length(args) >= 1) args[1] else "private_out"
POPFILE <- if (length(args) >= 2) args[2] else ""

files <- list.files(OUTDIR, pattern = "^private\\..*\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("private.<pop>.txt not found, please complete the filtering task first")
pop_of <- function(f) sub("^private\\.(.*)\\.txt$", "\\1", basename(f))
pops   <- pop_of(files)

# Private number per group (number of lines from file)
n_priv <- vapply(files, function(f)
  nrow(read.table(f, header = TRUE, sep = "\t")), integer(1))
names(n_priv) <- pops

# Number of candidates/number of samples (from stats file, if present)
sf  <- list.files(OUTDIR, pattern = "\\.stats\\.tsv$", full.names = TRUE)
if (length(sf) > 0) {
  st <- do.call(rbind, lapply(sf, function(f)
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
  st$safe <- gsub("[^A-Za-z0-9_.-]", "_", st$pop)
  summary <- data.frame(pop = st$pop, n_ind = st$n_ind,
                        n_candidate = st$n_candidate,
                        n_private = n_priv[st$safe],
                        stringsAsFactors = FALSE)
} else {
  summary <- data.frame(pop = pops, n_ind = NA_integer_,
                        n_candidate = NA_integer_,
                        n_private = n_priv[pops],
                        stringsAsFactors = FALSE)
}
write.table(summary, file.path(OUTDIR, "summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Cross-population deduplication of total private sites
all <- lapply(files, function(f) read.table(f, header = TRUE, sep = "\t",
                                            stringsAsFactors = FALSE))
all <- do.call(rbind, all)
all <- all[!duplicated(all$ID), ]
num_chr <- function(x) { n <- suppressWarnings(as.numeric(x)); ifelse(is.na(n), Inf, n) }
all <- all[order(num_chr(all$CHROM), all$POS), ]
write.table(all, file.path(OUTDIR, "private_multipop.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

png(file.path(OUTDIR, "summary.png"), width = 1600, height = 900, res = 150)
par(mar = c(10, 4.5, 2, 1))
b <- barplot(summary$n_private, names.arg = summary$pop,
             las = 2, col = "steelblue", border = NA,
             ylab = "Private SNPs", main = "Private SNPs per population")
text(b, summary$n_private, summary$n_private, pos = 3, cex = 0.7, col = "grey20")
dev.off()

cat(sprintf("Summary: %d populations, total number of private SNPs (duplicate removed): %d\n",
            nrow(summary), nrow(all)))
cat("Output: summary.tsv / private_multipop.tsv / summary.png\n")
