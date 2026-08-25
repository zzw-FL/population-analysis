#!/usr/bin/env Rscript
# =====================================================================
# p2_split_strat.R - Split the entire strat table into 1 file per group + 1 total group table
#
# Purpose: Each subsequent group task only reads 2 small tables, memory = O (number of sites) instead of O (site x population),
#       Each group can be safely run in parallel as an independent qsub task physically
#
# Usage: Rscript p2_split_strat.R <strat> <splits_dir>
# Output: splits/pop.<pop>.tsv per population (CHROM POS ID REF ALT ALT_CT OBS_CT)
#       splits/total.tsv Total population (CHROM POS ID TOT_ALT TOT_OBS)
# =====================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript p2_split_strat.R <strat> <splits_dir>")
STRAT  <- args[1]
SPLITS <- args[2]
dir.create(SPLITS, showWarnings = FALSE, recursive = TRUE)

has_dt <- requireNamespace("data.table", quietly = TRUE)
cat(sprintf("Reading %s ...\n", STRAT))
if (has_dt) {
  d <- data.table::fread(STRAT, check.names = FALSE, showProgress = FALSE)
  d <- as.data.frame(d, stringsAsFactors = FALSE)
} else {
  d <- read.table(STRAT, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, sep = "\t", comment.char = "")
}
names(d) <- sub("^#", "", names(d))

find_col <- function(aliases) {
  hit <- match(aliases, names(d)); hit <- hit[!is.na(hit)]
  if (length(hit) == 0) {
    stop(sprintf("Column [%s] not found. Actual column name: %s",
                 paste(aliases, collapse = "/"), paste(names(d), collapse = ", ")))
  }
  hit[1]
}
i_chrom <- find_col("CHROM");      i_pos <- find_col("POS")
i_id    <- find_col(c("ID", "SNP")); i_ref <- find_col("REF")
i_alt   <- find_col("ALT");        i_clu <- find_col(c("CLUSTER", "CLUST"))
i_act   <- find_col(c("ALT_CT", "A1_CT")); i_oct <- find_col(c("OBS_CT", "NCHROBS"))

chrom <- d[[i_chrom]]; pos <- d[[i_pos]]; idv <- as.character(d[[i_id]])
ref   <- d[[i_ref]];   alt <- d[[i_alt]]; clu <- as.character(d[[i_clu]])
act   <- as.numeric(d[[i_act]]); oct <- as.numeric(d[[i_oct]])

pops <- unique(clu)
safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
cat(sprintf("Split %d groups...\n", length(pops)))

# 1) One file per group
for (g in pops) {
  m  <- clu == g
  out <- data.frame(CHROM = chrom[m], POS = pos[m], ID = idv[m],
                    REF = ref[m], ALT = alt[m], ALT_CT = act[m], OBS_CT = oct[m],
                    stringsAsFactors = FALSE)
  f <- file.path(SPLITS, sprintf("pop.%s.tsv", safe(g)))
  if (has_dt) {
    data.table::fwrite(out, f, sep = "\t")
  } else {
    write.table(out, f, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  cat(sprintf("[%s] %d lines -> %s\n", g, nrow(out), basename(f)))
}

# 2) Total population table (aggregated by site)
uid  <- unique(idv)
idx  <- match(idv, uid)
tot_act <- rowsum(act, idx)[, 1]
tot_oct <- rowsum(oct, idx)[, 1]
first <- !duplicated(idv)
total <- data.frame(CHROM = chrom[first], POS = pos[first], ID = uid,
                    TOT_ALT = tot_act, TOT_OBS = tot_oct,
                    stringsAsFactors = FALSE)
ft <- file.path(SPLITS, "total.tsv")
if (has_dt) {
  data.table::fwrite(total, ft, sep = "\t")
} else {
  write.table(total, ft, sep = "\t", row.names = FALSE, quote = FALSE)
}
cat(sprintf("Total population table: %s (%d rows)\n", ft, nrow(total)))
