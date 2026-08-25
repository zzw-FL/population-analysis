#!/usr/bin/env Rscript
# =====================================================================
# Extract the CV errors in the ADMIXTURE log, summarize and plot them, and recommend the best K
# Usage: Rscript extract_cv.R [logdir] [pattern] [outprefix]
#   logdir log directory, default "." (current directory)
#   Pattern log file regular pattern (pattern of list.files), scans all files by default
#             By default, it is only recognized in subfolders such as admixture1-10
#             log.admixture.k*.*.out log (repeated numbers are automatically identified from the path);
#             Recurse the entire logdir only if there is no admixture subdirectory
#   outprefix output prefix, default "cv"
# Input file: ADMIXTURE run log/standard output, which contains
#           Lines like "CV error (K=5): 0.48371" (qsub's .o file or redirection log)
# Output:
#   cv_raw.txt Each line: K repeat sequence number CV value log file
#   cv_summary.txt  K  n  mean  sd  min  max
#   cv_plot.pdf/png K-CV line chart (point=mean, error bar=min~max, mark the best K)
# Best K Recommendations:
#   bestK_min = K with the smallest CV mean (ADMIXTURE official standard)
#   bestK_parsimony = smallest K that is < 1% different from the smallest CV (more parsimonious, easier to interpret)
# Dependencies: base R only
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
logdir    <- if (length(args) >= 1) args[1] else "."
pattern   <- if (length(args) >= 2) args[2] else NULL
outprefix <- if (length(args) >= 3) args[3] else "cv"

# Default: only log.admixture.k*.*.out log files are recognized in the admixture{N} subdirectory
# (Fast and will not misread other files); Recurse the entire logdir when there is no such subdirectory;
# If pattern is given, match pattern
LOG_PAT <- "^log\\.admixture\\.k.*\\.out$"
if (is.null(pattern)) {
  subdirs <- list.dirs(logdir, recursive = FALSE, full.names = TRUE)
  subdirs <- subdirs[grepl("admixture[0-9]+$", subdirs)]
  if (length(subdirs) > 0) {
    files <- unlist(lapply(subdirs, list.files,
                           pattern = LOG_PAT, full.names = TRUE), use.names = FALSE)
  } else {
    files <- list.files(logdir, pattern = LOG_PAT, full.names = TRUE, recursive = TRUE)
  }
} else {
  files <- list.files(logdir, pattern = pattern, full.names = TRUE, recursive = TRUE)
}
# Exclude its own output files to avoid using the results as input during the second run.
files <- files[!grepl("\\.(pdf|png|txt)$", files, ignore.case = TRUE)]
# Read only ordinary small files: skip directories and large files >10MB (such as .bed/.phy),
# Avoid readLines freezing or running out of memory (ADMIXTURE logs are usually only a few KB)
info <- file.info(files)
files <- files[!is.na(info$size) & !info$isdir & info$size <= 10 * 1024^2]
cat("Scanning", length(files), "log files...\n")

cv_list <- list()
for (f in files) {
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) next
  m <- regexec("CV error \\(K=([0-9]+)\\): ([0-9.]+)", lines)
  hits <- regmatches(lines, m)
  for (h in hits) {
    if (length(h) == 3) {
      cv_list[[length(cv_list) + 1]] <-
        data.frame(K = as.integer(h[2]), CV = as.numeric(h[3]),
                   file = sub(paste0("^", logdir, "/+"), "", f),
                   stringsAsFactors = FALSE)
    }
  }
}
if (length(cv_list) == 0)
  stop("'CV error (K=..): ..' line not found in any logs, please check logdir/pattern")

raw <- do.call(rbind, cv_list)
# Duplicate numbers are first extracted from the admixture{N} folder in the path; if not, they are numbered in the order in which they appear within each K.
rep_p <- sapply(raw$file, function(x) {
  hh <- regmatches(x, regexec("admixture([0-9]+)", x))[[1]]
  if (length(hh) == 2) as.integer(hh[2]) else NA_integer_
})
if (all(is.na(rep_p))) {
  raw$rep <- ave(seq_len(nrow(raw)), raw$K, FUN = seq_along)
} else {
  raw$rep <- rep_p
}

# ---------- Summary ----------
summ <- aggregate(CV ~ K, raw, function(x)
  c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)))
summ <- do.call(data.frame, summ)
names(summ) <- c("K", "n", "mean", "sd", "min", "max")
summ <- summ[order(summ$K), ]

write.table(raw[, c("K", "rep", "CV", "file")], sprintf("%s_raw.txt", outprefix),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(summ, sprintf("%s_summary.txt", outprefix),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----------Best K ----------
cv_min <- min(summ$mean)
best_min <- summ$K[which.min(summ$mean)]
tol <- 0.01 * cv_min
best_pars <- min(summ$K[summ$mean <= cv_min + tol])
top3 <- summ[order(summ$mean), ][1:min(3, nrow(summ)), c("K", "mean")]

cat("========== CV error summary (repeated mean for each K) ==========\n")
print(summ, row.names = FALSE)
cat(sprintf("\nBest K (minimum CV): K = %d (CV = %.4f)\n", best_min, cv_min))
cat(sprintf("Best K (parsimonious, 1%% tolerance): K = %d (CV = %.4f)\n",
            best_pars, summ$mean[summ$K == best_pars]))
cat("Top 3 K's with lowest CV:\n")
print(top3, row.names = FALSE)
cat(sprintf("\nRecommendation: K = %d is preferred; if the CV difference between adjacent K's is <1%% and the biological explanation is more reasonable, K = %d can be used instead\n",
            best_min, best_pars))

# ---------- Drawing ----------
draw_cv <- function() {
  par(mar = c(4.5, 4.5, 3, 1))
  plot(summ$K, summ$mean, type = "b", pch = 19, col = "steelblue", lwd = 2,
       xlab = "K", ylab = "CV error",
       ylim = range(c(summ$min, summ$max)),
       main = "ADMIXTURE cross-validation error")
  arrows(summ$K, summ$min, summ$K, summ$max,
         angle = 90, code = 3, length = 0.04, col = "grey50")
  abline(v = best_min, col = "red", lty = 2)
  points(best_min, cv_min, col = "red", pch = 19, cex = 1.6)
  text(best_min, max(summ$max), labels = sprintf("best K = %d", best_min),
       col = "red", pos = 3, xpd = NA)
}
pdf(sprintf("%s_plot.pdf", outprefix), width = 7, height = 5)
draw_cv()
dev.off()
png(sprintf("%s_plot.png", outprefix), width = 1400, height = 1000, res = 150)
draw_cv()
dev.off()

cat(sprintf("\nOutput: %s_raw.txt %s_summary.txt %s_plot.pdf %s_plot.png\n",
            outprefix, outprefix, outprefix, outprefix))
