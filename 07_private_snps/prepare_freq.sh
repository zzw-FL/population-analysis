#!/bin/bash
# =====================================================================
# p1_prepare_freq.sh - Single full scan, generate site x population allele count table
#
# Input: large group VCF (or converted BED), popfile (sample ID<TAB> group)
# Output: ${OUT}.acount.strat[.gz] (plink2 --freq counts --within results)
#
# usage:
#   VCF=cycas_whole_genome.vcf.gz bash p1_prepare_freq.sh # Scan VCF directly
#   BED=cycas.maffiltered.4dtv.snp bash p1_prepare_freq.sh # Scan BED (faster)
# It is recommended to submit the cluster: qsub -cwd -pe smp 4 p1_prepare_freq.sh
#
# Configurable: PLINK2 / VCF / BED / POPFILE / OUT / THREADS / SKIP
#   PLINK2 plink2 path, such as /TOOLS/USER/software/plink2/plink2
#   POPFILE two columns: sample ID<TAB> population, no header (if there is a header, set SKIP=1)
# =====================================================================
set -euo pipefail

PLINK2="/TOOLS/plink2"
VCF="/DATA/01.vcf/03.maf.merge/SAMPLE.10chr.maf.snp.vcf.gz"
BED="/DATA/01.vcf/03.maf.merge/SAMPLE.10chr.maf.snp"
POPFILE="sample.list"
OUT="uniq_SNP"
THREADS="4"
SKIP="${SKIP:-0}"

command -v "$PLINK2" >/dev/null 2>&1 || { echo "Error: plink2 not found, please set PLINK2=..." >&2; exit 1; }
[ -f "$POPFILE" ] || { echo "Error: Missing $POPFILE" >&2; exit 1; }
if [ -n "$VCF" ]; then
  [ -f "$VCF" ] || { echo "Error: Missing $VCF" >&2; exit 1; }
else
  [ -f "$BED.bed" ] || { echo "Error: Missing $BED.bed (or set VCF=xxx.vcf.gz)" >&2; exit 1; }
fi

# 1) popfile -> plink --within three-column file (FID IID CLUSTER, FID=IID)
if [ "$SKIP" = "1" ]; then
  tail -n +2 "$POPFILE" > "${POPFILE}.nohdr"; SRC="${POPFILE}.nohdr"
else
  SRC="$POPFILE"
fi
awk -F '\t' 'NF==2 {print $1"\t"$1"\t"$2} NF>=3 {print $1"\t"$2"\t"$3}' "$SRC" > "${POPFILE}.within"
[ "$SKIP" = "1" ] && rm -f "${POPFILE}.nohdr"

# 2) Single scan: Output ALT copy number and observed allele number for each variant and each population
if [ -n "$VCF" ]; then
  echo "Single scan from VCF: $VCF"
  "$PLINK2" --vcf "$VCF" --double-id --chr-set 33 --set-missing-var-ids @:# \
    --within "${POPFILE}.within" --freq counts \
    --threads "$THREADS" --out "$OUT"
else
  echo "Single scan from BED: $BED"
  "$PLINK2" --bfile "$BED" --chr-set 33 --within "${POPFILE}.within" \
    --freq counts --threads "$THREADS" --out "$OUT"
fi

# 3) Check the output
STRAT="${OUT}.acount.strat.gz"
[ -f "$STRAT" ] || STRAT="${OUT}.acount.strat"
[ -f "$STRAT" ] || { echo "Error: ${OUT}.acount.strat.gz not found" >&2; exit 1; }
echo "Done: $STRAT"
echo "Next step: Rscript p2_split_strat.R $STRAT splits"
