#!/bin/bash
# =====================================================================
# p3_make_jobs.sh - generate filter script + batch submission script for each group to submit clusters
#
# Usage: bash p3_make_jobs.sh
# Generate: extract_private.<pop>.sh (one for each group, can be directly qsub)
#       submit_all.sh (one-click submission of all groups, physical parallelism)
# Configurable: POPFILE / SPLITS / OUTDIR / MAC / MIN_GENO / OTHERS_AF_MAX / PE
# =====================================================================
set -euo pipefail

POPFILE="${POPFILE:-sample.pop.txt}"
SPLITS="${SPLITS:-splits}"
OUTDIR="${OUTDIR:-private_out}"
MAC="${MAC:-2}"; MIN_GENO="${MIN_GENO:-0.8}"; OTHERS_AF_MAX="${OTHERS_AF_MAX:-0}"
PE="${PE:-smp 2}"
ROOT="$(cd "$(dirname "$0")" && pwd)"

[ -f "$POPFILE" ] || { echo "Error: Missing $POPFILE" >&2; exit 1; }
[ -f "$SPLITS/total.tsv" ] || { echo "Error: Missing $SPLITS/total.tsv, run p2 first" >&2; exit 1; }

n=0
for pop in $(cut -f2 "$POPFILE" | sort -u); do
  safe=$(echo "$pop" | tr ' /()' '____')
  [ -f "$SPLITS/pop.${safe}.tsv" ] || { echo "Warning: Missing $SPLITS/pop.${safe}.tsv, skipping $pop" >&2; continue; }
  cat > "extract_private.${safe}.sh" <<EOF
#!/bin/bash
Rscript $ROOT/p4_filter_pop.R $SPLITS/pop.${safe}.tsv $SPLITS/total.tsv $POPFILE $OUTDIR "$pop" $MAC $MIN_GENO $OTHERS_AF_MAX
EOF
  chmod +x "extract_private.${safe}.sh"
  n=$((n+1))
done

cat > submit_all.sh <<EOF
#!/bin/bash
# Submit all group tasks with one click (physical parallelism); run p5 summary after all are completed
for s in extract_private.*.sh; do
  qsub -cwd -pe $PE -N "priv_\${s%.sh}" "\$s"
done
EOF
chmod +x submit_all.sh

echo "Generate $n group task scripts; submit: bash submit_all.sh"
echo "After all is completed: Rscript p5_summarize.R $OUTDIR $POPFILE"
