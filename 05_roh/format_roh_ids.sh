#!/bin/bash
set -e

# ==================== 1. Configuration area ====================
VCF_DIR="/DATA/01.vcf/05.population_maf"
VCF_SUFFIX=".list.maf.snp.vcf.gz"
PLINK="/TOOLS/USER/software/plink/plink"

# Explicitly specify 5 species that previously failed
FAILED_POPS=(
    "C.longiconifera"
    "C.micronesica"
    "C.multipinnata"
    "C.sexseminifera"
    "C.taiwaniana"
)
# ===================================================

echo "========== Phase 1: Preparing for serial reruns of 5 failed species =========="
REPORT="Rerun_ROH_Summary.txt"
# Clear or create a new results file
> $REPORT

for POP_NAME in "${FAILED_POPS[@]}"; do
    VCF_FILE="${VCF_DIR}/${POP_NAME}${VCF_SUFFIX}"

    if [ ! -f "$VCF_FILE" ]; then
        echo "⚠️ File not found: $VCF_FILE, skipping..."
        continue
    fi

    echo "▶️ Concentrating all resources for calculation: $POP_NAME ..."

    MEM_MB=50000

    # 1. Convert format (added --set-missing-var-ids to fix ID pitfall prevention)
    $PLINK --vcf "$VCF_FILE" \
           --double-id --allow-extra-chr --chr-set 34 \
           --set-missing-var-ids @:# \
           --memory $MEM_MB \
           --make-bed \
           --out "tmp_${POP_NAME}" > /dev/null 2>&1

    # 2. Read the binary file and run ROH
    $PLINK --bfile "tmp_${POP_NAME}" \
           --allow-extra-chr --chr-set 34 \
           --memory $MEM_MB \
           --homozyg \
           --homozyg-kb 100 \
           --homozyg-snp 50 \
           --homozyg-window-snp 50 \
           --homozyg-window-het 3 \
           --homozyg-window-missing 50 \
           --out "roh_${POP_NAME}" > /dev/null 2>&1

    # 3. Check and extract results
    if [ -f "roh_${POP_NAME}.hom.indiv" ]; then
        SUMMARY=$(awk '
        NR>1 {
            nseg += $4;
            total_kb += $5;
            count++;
        }
        END {
            if(count>0)
                printf "%d\t%.2f\t%.2f\n", count, nseg/count, (total_kb/1000)/count;
            else
                print "0\t0.00\t0.00";
        }' "roh_${POP_NAME}.hom.indiv")

        echo -e "${POP_NAME}\t${SUMMARY}" >> $REPORT
    else
        echo -e "${POP_NAME}\t0\t0.00\t0.00\tError" >> $REPORT
    fi

    # 4. Clean the battlefield
    rm -f "tmp_${POP_NAME}.bed" "tmp_${POP_NAME}.bim" "tmp_${POP_NAME}.fam" "tmp_${POP_NAME}.log" "tmp_${POP_NAME}.nosex"

    echo "✅ $POP_NAME done!"
    echo "-----------------------------------"

done

echo -e "\n========== 🎉 The rerun results are as follows =========="
echo -e "Species/Population\tN_Samples\tMean_ROH_Count\tMean_Total_Length(Mb)"
cat $REPORT | column -t
echo "--------------------------------------------------------"
echo "You can directly copy and paste these 5 rows of results into the original large table!"
