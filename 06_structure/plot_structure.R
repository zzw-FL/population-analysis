#!/usr/bin/env Rscript
# =====================================================================
# CLUMPP results plotting: Structure style ancestor component bar chart (relies on base R only)
# Input: admixture_k{K}.outfile CLUMPP average Q matrix, supports three formats:
#         1) Pure Q matrix (K columns)
#         2) Individual name in the first column + K column Q (K+1 column)
#         3) CLUMPP native output: serial number individual number (0) cluster number: Q1...QK (K+5 columns)
#       Optional admixture_k{K}.indfile (individual name, supports CLUMPP "#1 name" format),
#            popfile.txt (individual name group),
#            sample.order.list (ML tree sample order, one sample name per line)
# Sorting: The default (auto) priority is to sort by sample.order.list (ML tree order),
#       When there is no such file: there is popfile grouping by population, otherwise clustering by Q value
# Color: Nature (NPG) series color matching: The first 10 clusters use Nature 10 colors, and the remaining parts use
#       CIELAB perceptual distance maximization for color selection from the Nature hue family extension; automatic clustering across K
#       Matching: The color of each old cluster is only passed to the branch with the "largest color block" among its descendant clusters.
#       (The samples whose dominant cluster is j in the previous K, and the confusion number whose dominant cluster is i in the current K
#       C[i,j] is the largest cluster), when splitting, large blocks keep the old color, and new small blocks are assigned new colors;
#       The order of color numbers is to keep the old color on the left and the new color after it according to size to reduce visual jitter as much as possible.
# Output: structure_K{K}.png Each K single image
#       structure_allK.pdf K range splicing large picture (each row is sorted by its own K)
#       structure_comb_K2-K4-K8.pdf/png specifies K combination mosaic diagram (combineK parameter)
# Usage: Rscript plot_structure.R [Kmin] [Kmax] [popfile] [sortmode] [refK] [orderfile] [combineK]
# Example: Rscript plot_structure.R 2 35 # auto: orderfile > pop > q
#       Rscript plot_structure.R 2 35 popfile.txt # Same as above
#       Rscript plot_structure.R 2 35 popfile.txt pop # Force grouping
#       Rscript plot_structure.R 2 35 popfile.txt q 8 # Sort fixed panels with K=8
#       Rscript plot_structure.R 2 35 - order - tree.order.list # Specify the order file
#       Rscript plot_structure.R 2 35 popfile.txt auto - - "2,4,6,8,10-24" # Additional output specifies K combination plot
# sortmode: auto(default) | order | q | dom | pop | name
#   order = sort by sample.order.list (unlisted samples are at the end)
#   q = Ranking based on Q matrix hierarchical clustering (similar genetic backgrounds cluster together)
#   dom = sort by dominant cluster + dominant proportion in descending order
#   pop = group by group, sort by dominant cluster within the group
#   name = Sort by individual name (pure numbers by numerical value)
# refK: Optional, fix all panels with the sorting results of the specified K (multi-line alignment), only valid for q/dom
# orderfile: optional, default sample.order.list; pass "-" to disable
# combineK: optional, comma separated list of K, supports range (such as "2,4,6,8,10-24");
#           Additional output is the combined mosaic image in the same style as allK structure_comb_*.pdf/png
# =====================================================================

args       <- commandArgs(trailingOnly = TRUE)
Kmin       <- if (length(args) >= 1) as.integer(args[1]) else 2
Kmax       <- if (length(args) >= 2) as.integer(args[2]) else 35
popfile    <- if (length(args) >= 3) args[3] else "popfile.txt"
sortmode   <- if (length(args) >= 4) tolower(args[4]) else "auto"
refK       <- if (length(args) >= 5 && nzchar(args[5]) && args[5] != "-") suppressWarnings(as.integer(args[5])) else NA
orderfile  <- if (length(args) >= 6) args[6] else "sample.order.list"
combineK   <- if (length(args) >= 7) args[7] else ""

# ---------- Analysis combineK: "2,4,6,8,10-24" -> c(2,4,6,8,10:24) ----------
parse_combine <- function(s) {
  s <- gsub("[[:space:]]", "", s)
  if (!nzchar(s)) return(integer(0))
  out <- integer(0)
  for (tok in unlist(strsplit(s, ",", fixed = TRUE))) {
    if (grepl("^[0-9]+-[0-9]+$", tok)) {
      ab <- as.integer(strsplit(tok, "-", fixed = TRUE)[[1]])
      out <- c(out, seq(ab[1], ab[2]))
    } else if (grepl("^[0-9]+$", tok)) {
      out <- c(out, as.integer(tok))
    } else {
      cat("Warning: Unable to parse K combination '", tok, "', ignored\n", sep = "")
    }
  }
  sort(unique(out))
}
comb_ks <- parse_combine(combineK)


# ---------- Color: Nature (NPG) series color matching, supports dozens of clusters ----------
# Skeleton: Nature Publishing Group 10 colors (ggsci::pal_npg):
#   #E64B35 #4DBBD5 #00A087 #3C5488 #F39B7F
#   #8491B4 #91D1C2 #DC0000 #7E6148 #B09C85
# K<=10 uses these 10 colors directly in order (consistent with the Nature style map);
# When K>10, the new cluster is maximized according to CIELAB perceptual distance, from the Nature hue family
# (10 anchor point hue x 3 levels of brightness x 2 levels of chroma x hue fine-tuning) Extended color selection,
# It not only maintains the Nature texture but also ensures that they can be separated from each other.
nature_base <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                 "#7E6148", "#B09C85")

col2lab <- function(cols) {
  rgb <- t(col2rgb(cols)) / 255
  convertColor(rgb, from = "sRGB", to = "Lab")
}
# Lab -> HCL polar coordinates (Hue, Chroma, Lightness), extract the hue anchor point of the Nature base color
lab2hcl <- function(cols) {
  lab <- col2lab(cols)
  h <- atan2(lab[, 3], lab[, 2]) * 180 / pi
  h[h < 0] <- h[h < 0] + 360
  cbind(h = h, c = sqrt(lab[, 2]^2 + lab[, 3]^2), l = lab[, 1])
}
base_hcl <- lab2hcl(nature_base)
# Nature Hue Family Candidate Palette: 10 anchor hue x 3 hue offset x 3 lightness x 2 chroma = 180 colors
hue_offsets <- c(-18, 0, 18)
l_lev <- c(72, 55, 38)
c_lev <- c(95, 68)
cand_grid <- expand.grid(h = base_hcl[, "h"], off = hue_offsets,
                         l = l_lev, c = c_lev)
palette_cand_hex <- hcl(cand_grid$h + cand_grid$off,
                        c = cand_grid$c, l = cand_grid$l)
palette_cand_lab <- col2lab(palette_cand_hex)
pick_new_colors <- function(used_lab, n) {
  cand <- seq_along(palette_cand_hex)
  out <- integer(0)
  for (k in seq_len(n)) {
    ref <- rbind(used_lab, palette_cand_lab[out, , drop = FALSE])
    d <- apply(palette_cand_lab[cand, , drop = FALSE], 1, function(cc)
      min(sqrt(rowSums((ref - rep(cc, each = nrow(ref)))^2))))
    i <- cand[which.max(d)]
    out <- c(out, i)
    cand <- setdiff(cand, i)
  }
  palette_cand_hex[out]
}
# Color matching of first K clusters: K<=10 uses Nature 10 colors; K>10 extends beyond Nature base color
q_colors <- function(n) {
  if (n <= length(nature_base)) return(nature_base[seq_len(n)])
  c(nature_base, pick_new_colors(col2lab(nature_base), n - length(nature_base)))
}

# ---------- Sorting function ----------
# 1) According to Q similarity: hierarchical clustering, individuals with similar genetic backgrounds are arranged together
order_q <- function(q) {
  if (nrow(q) < 3) return(seq_len(nrow(q)))
  hclust(dist(q), method = "average")$order
}
# 2) According to dominant cluster + dominant proportion: first group by the cluster to which the largest Q belongs, and then sort in descending order according to proportion within the cluster
order_dom <- function(q) {
  if (nrow(q) < 2) return(seq_len(nrow(q)))
  dom <- max.col(q, ties.method = "first")
  order(dom, -apply(q, 1, max))
}
# 3) Sort individual names: pure numbers are sorted by numerical value, otherwise they are sorted in dictionary order
num_order <- function(x) {
  nums <- suppressWarnings(as.numeric(x))
  if (all(!is.na(nums))) order(nums) else order(x)
}
# 4) File sorting by sample order: samples appearing in ord_ids are in file order, and the rest are sorted at the end
order_list <- function(ids, ord_ids) {
  ord <- match(ord_ids, ids)
  ord <- ord[!is.na(ord)]
  if (length(ord) == 0) return(NULL)
  c(ord, setdiff(seq_along(ids), ord))
}
# 5) Cross-K color inheritance: Use the confusion matrix of the sample-dominant cluster to match adjacent K clusters,
#    Return perm[i] = the color number that the i-th cluster of current K should use
#    (The length is always equal to Kc, no columns are lost; the matching cluster uses the color number of the previous K, and the new cluster uses the new number of Kp+1..Kc)
#    Kc/Kp is the real number of columns of the current/previous K, not the number of unique(dom):
#    A cluster must retain its Q column even if it does not have a dominant sample, otherwise white blocks will appear when the row sum <1
#    Matching rules: The color of each old cluster j is inherited by its descendant cluster with the "largest color block", that is, in the previous K
#    The sample whose dominant cluster is j has the largest confusion number C[i,j] of the dominant cluster i in the current K;
#    In the event of a tie, the leftmost cluster is taken (the position is more stable); when the old cluster is split, the large block retains the old color.
#    Newly appearing small blocks are assigned new color numbers, and the new clusters are arranged after the old colors in descending order of size.
match_clusters <- function(dom_prev, dom_curr, Kc, Kp) {
  C <- matrix(0L, nrow = Kc, ncol = Kp)
  for (i in seq_len(Kc)) {
    for (j in seq_len(Kp)) {
      C[i, j] <- sum(dom_curr == i & dom_prev == j)
    }
  }
  # Panel position of the leftmost sample of each current cluster (only used for position stabilization in the case of a tie)
  leftmost <- vapply(seq_len(Kc), function(i) {
    w <- which(dom_curr == i)
    if (length(w) > 0) min(w) else Inf
  }, numeric(1))

  perm <- integer(Kc)
  used_i <- logical(Kc)
  # Process the previous K clusters one by one: select the one with the largest confusion number among its descendant clusters (C[,j]>0 and not allocated),
  # That is, the largest color block inherits the original color number; in a tie, the leftmost cluster is selected.
  for (j in seq_len(Kp)) {
    cand <- which(C[, j] > 0L & !used_i)
    if (length(cand) == 0) next
    i <- cand[order(-C[cand, j], leftmost[cand])[1]]
    perm[i] <- j
    used_i[i] <- TRUE
  }
  # Current cluster on no match (new cluster/no dominant sample cluster): assign new color numbers in descending order of size
  rem <- which(!used_i)
  if (length(rem) > 0) {
    sizes <- tabulate(dom_curr, Kc)
    rem <- rem[order(sizes[rem], decreasing = TRUE)]
    perm[rem] <- Kp + seq_along(rem)
  }
  perm
}

# ---------- Reading sample order (ML tree order, optional) ----------
ord_ids <- NULL
if (orderfile != "-" && file.exists(orderfile)) {
  ord_ids <- trimws(readLines(orderfile))
  ord_ids <- ord_ids[nzchar(ord_ids)]
}
has_order <- length(ord_ids) > 0

# ---------- Read the individual name (any indfile of K is enough, the content should be consistent) ----------
ids <- NULL
for (K in Kmin:Kmax) {
  f <- sprintf("admixture_k%d.indfile", K)
  if (file.exists(f)) {
    ids <- trimws(readLines(f))
    ids <- sub("^# [0-9]+[[:space:]]*", "", ids) # Compatible with CLUMPP "#1 name" format
    break
  }
}

# ---------- Read group information (optional) ----------
# popfile two columns: individual name <TAB> group name; or single column (group names are given in order of individuals)
pop <- NULL
if (file.exists(popfile)) {
  p <- read.table(popfile, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(p) == 1 && !is.null(ids) && nrow(p) == length(ids)) {
    pop <- setNames(p[[1]], ids)
  } else if (ncol(p) >= 2) {
    pop <- setNames(p[[2]], p[[1]])
  }
}
has_pop <- !is.null(pop)
if (sortmode == "auto") {
  if (has_order) sortmode <- "order"
  else if (has_pop) sortmode <- "pop"
  else sortmode <- "q"
}
if (sortmode == "order" && !has_order) {
  cat("Warning: No available sample.order.list found, fallback to q sort\n")
  sortmode <- "q"
}
if (sortmode == "pop" && !has_pop) {
  cat("Warning: No available popfile found, fallback to q sort\n")
  sortmode <- "q"
}

# ---------- Drawing function: one K, one panel ----------
draw_panel <- function(q, cols, seps = NULL, xaxs = "r", yaxs = "r") {
  n <- nrow(q); K <- ncol(q)
  plot.new()
  plot.window(xlim = c(0, n), ylim = c(0, 1), xaxs = xaxs, yaxs = yaxs)
  for (i in seq_len(n)) {
    cum <- 0
    for (j in seq_len(K)) {
      h <- q[i, j]
      if (h > 0) {
        rect(i - 1, cum, i, cum + h, col = cols[j], border = NA)
        cum <- cum + h
      }
    }
  }
  if (!is.null(seps)) {
    abline(v = seps, lty = 2, col = "grey40", lwd = 0.5)
  }
}

done <- c()
ref_order <- NULL   # Sorted results of refK, used to pin all panels
dom_prev  <- NULL   # Sample-dominant clusters for the last K (for color inheritance matching)
K_prev    <- NULL   # Number of columns for the previous K (for color inheritance matching)
cols_prev <- NULL   # The actual color of the previous K (matching clusters are strictly used to avoid color drift)

for (K in Kmin:Kmax) {
  f <- sprintf("admixture_k%d.outfile", K)
  if (!file.exists(f)) { cat("Skip (no file):", f, "\n"); next }

  q <- read.table(f, header = FALSE, stringsAsFactors = FALSE)
  nc <- ncol(q)
  ids2 <- NULL
  if (nc == K) {
    ids2 <- ids                          # Pure Q matrix
  } else if (nc == K + 1) {
    ids2 <- as.character(q[[1]])         # First column individual name
    q <- q[, -1, drop = FALSE]
  } else if (nc == K + 5) {
    # CLUMPP native output: serial number individual number (0) cluster number: Q1...QK
    ids2 <- as.character(q[[2]])         # Column 2 individual number
    q <- q[, (nc - K + 1):nc, drop = FALSE]
  } else {
    cat("skip (actual number of columns", nc, "!= K", K, "):", f, "\n")
    next
  }
  q <- as.matrix(q)
  storage.mode(q) <- "numeric"

  # Individual identification: Priority is given to the indfile name (the same as the line order of outfile), otherwise the outfile name is used.
  if (!is.null(ids) && length(ids) == nrow(q)) {
    rownames(q) <- ids
  } else if (!is.null(ids2) && length(ids2) == nrow(q)) {
    rownames(q) <- ids2
  }

  # ---------- Sample sorting ----------
  ord <- seq_len(nrow(q))
  if (sortmode == "name") {
    if (!is.null(rownames(q))) ord <- num_order(rownames(q))
  } else if (sortmode == "dom") {
    ord <- order_dom(q)
  } else if (sortmode == "order") {
    # When the indfile sample name is inconsistent with orderfile, try to use outfile's own individual name instead.
    if (!is.null(rownames(q)) && length(intersect(rownames(q), ord_ids)) == 0 &&
        !is.null(ids2) && length(intersect(ids2, ord_ids)) > 0) {
      rownames(q) <- ids2
      cat("Tip: indfile sample name is inconsistent with orderfile, outfile individual name has been used instead\n")
    }
    if (is.null(rownames(q))) {
      cat("Warning: Unable to obtain sample name, fallback to q sort\n")
      ord <- order_q(q)
    } else {
      ord <- order_list(rownames(q), ord_ids)
      if (is.null(ord)) {
        cat("Warning: sample name cannot match sample.order.list, fallback to q order\n")
        ord <- order_q(q)
      } else {
        n_miss <- nrow(q) - length(intersect(rownames(q), ord_ids))
        if (n_miss > 0) cat("Warning: orderfile not included", n_miss, "samples, sorted to the end\n")
        cat("  K=", K, ": Sort by sample.order.list\n", sep = "")
      }
    }
  } else if (sortmode == "pop") {
    if (has_pop && all(rownames(q) %in% names(pop))) {
      pv <- as.character(pop[rownames(q)])
      ord <- order(pv, order_dom(q))     # The group comes first, and the group is divided according to the dominant cluster + proportion
    } else {
      cat("Warning: Individual name cannot match popfile, fallback to q sort\n")
      ord <- order_q(q)
    }
  } else {                               # q
    ord <- order_q(q)
  }

  # refK: fixed ordering; reuse only when sample order can be matched (q/dom mode)
  if (!is.na(refK) && K == refK) {
    ref_order <- ord
  } else if (!is.na(refK) && sortmode %in% c("q", "dom") && !is.null(ref_order)) {
    ord <- ref_order
  }
  q <- q[ord, , drop = FALSE]

  # ---------- Cross-K color inheritance: rearrange Q columns after matching clusters ----------
  # The color of each old cluster is inherited by its descendant cluster with the "largest color block" (when splitting, large blocks retain the color and small blocks change to new colors),
  # Clusters with no dominant sample or no matches are always assigned a new color; the new color is from the Nature hue family
  # Maximize color selection according to CIELAB distance in candidate board
  dom_curr <- max.col(q, ties.method = "first")
  if (!is.null(dom_prev) && length(dom_prev) == nrow(q)) {
    perm <- match_clusters(dom_prev, dom_curr, K, K_prev)
    perm_s <- perm[order(perm)]        # Color number corresponding to column i after rearrangement
    q <- q[, order(perm), drop = FALSE]
    cols <- rep(NA_character_, K)
    for (i in seq_len(K)) {
      if (perm_s[i] <= K_prev) cols[i] <- cols_prev[perm_s[i]]   # The cluster on the left strictly uses the original color
    }
    new_n <- sum(is.na(cols))
    if (new_n > 0) {
      used_lab <- col2lab(cols[!is.na(cols)])
      cols[is.na(cols)] <- pick_new_colors(used_lab, new_n)  # New clusters with new colors
    }
  } else {
    sizes <- tabulate(dom_curr, K)
    q <- q[, order(sizes, decreasing = TRUE), drop = FALSE]  # First K: Big cluster first
    cols <- q_colors(K)
  }
  dom_prev <- max.col(q, ties.method = "first")
  K_prev <- K
  cols_prev <- cols

  pop_ord <- NULL
  seps <- NULL
  if (has_pop && all(rownames(q) %in% names(pop))) {
    pop_ord <- as.character(pop[rownames(q)])
    seps <- cumsum(rle(pop_ord)$lengths)
  }

  n <- nrow(q)

  # ---Single K chart ---
  png(sprintf("structure_K%d.png", K), width = 2800, height = 900, res = 300)
  par(mar = c(7, 4, 3, 1))
  draw_panel(q, cols, seps)
  if (!is.null(pop_ord)) {
    r <- rle(pop_ord)
    mids <- c(0, cumsum(r$lengths)[-length(r$lengths)]) + r$lengths / 2
    axis(1, at = mids, labels = r$values, las = 2, cex.axis = 0.6, tick = FALSE)
  }
  mtext(sprintf("K = %d", K), side = 3, adj = 0, line = 0.8, cex = 1.3)
  mtext("Individuals", side = 1, line = 5.5, cex = 0.9)
  mtext("Ancestry proportion", side = 2, line = 2.5, cex = 0.9)
  dev.off()
  cat("Output: structure_K", K, ".png (Sort by:", sortmode, ")\n", sep = "")

  done <- c(done, K)
  assign(sprintf("q%d", K), q)
  assign(sprintf("cols%d", K), cols)
  assign(sprintf("seps%d", K), seps)
}

# ---------- K range splicing large picture (shared drawing board function) ----------
draw_stack <- function(ks) {
  par(mfrow = c(length(ks), 1), mar = c(0, 5.5, 0, 1), oma = c(1.5, 0, 1, 0), xpd = NA)
  for (K in ks) {
    q  <- get(sprintf("q%d", K))
    se <- get(sprintf("seps%d", K))
    draw_panel(q, get(sprintf("cols%d", K)), se, xaxs = "i", yaxs = "i")
    mtext(sprintf("K = %d", K), side = 2, line = 1, las = 1, cex = 1.2)
  }
}
if (length(done) > 0) {
  pdf("structure_allK.pdf", width = 12, height = 0.45 * length(done) + 1)
  draw_stack(done)
  dev.off()
  cat("Output: structure_allK.pdf\n")
}

# ---------- Specify K combination mosaic (combineK, same as allK) ----------
if (length(comb_ks) > 0) {
  use <- comb_ks[comb_ks %in% done]
  if (length(use) > 0) {
    miss <- setdiff(comb_ks, done)
    if (length(miss) > 0)
      cat("Warning: No result file specified for K, skipped:", paste(miss, collapse = ","), "\n")
    stem <- paste("structure_comb", paste(use, collapse = "-"), sep = "_")
    pdf(sprintf("%s.pdf", stem), width = 12, height = 0.45 * length(use) + 1)
    draw_stack(use)
    dev.off()
    png(sprintf("%s.png", stem), width = 2400, height = (0.45 * length(use) + 1) * 200, res = 200)
    draw_stack(use)
    dev.off()
    cat("Output:", stem, ".pdf / ", stem, ".png\n", sep = "")
  } else {
    cat("Warning: There is no result file for the specified K combinations, and no combination graph is output\n")
  }
}
