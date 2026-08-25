#!/bin/bash
set -e

# ==================== 1. Configuration area ====================
# [Extremely important]: You must use pure MAF filter files that have not been LD thinned!
VCF_DIR="/DATA/01.vcf/05.population_maf"
VCF_SUFFIX=".list.maf.snp.vcf.gz"

PLINK="/TOOLS/USER/software/plink/plink"  # If you do not configure environment variables, write /absolute path/plink

# Maximum number of files allowed to be calculated simultaneously
# (Note: The script has limited each task to 50GB. If MAX_JOBS=6, a total of more than 300GB of memory must be applied for when submitting)
MAX_JOBS=6
# ===================================================

echo "========== Stage 1: Scan MAF VCF files =========="
VCF_FILES=(${VCF_DIR}/*${VCF_SUFFIX})
if [ ${#VCF_FILES[@]} -eq 0 ] || [ ! -e "${VCF_FILES[0]}" ]; then
    echo "❌ Error: No file ending with $VCF_SUFFIX found in directory!"
    exit 1
fi

FILE_COUNT=${#VCF_FILES[@]}
echo "✅ $FILE_COUNT VCF files successfully detected for ROH scan."

echo -e "\n========== Phase 2: Enable concurrent computing ROH =========="
REPORT="Final_Species_ROH_Summary.txt"
# The header of the output report: species name, number of samples, average number of ROH, average total length of ROH (Mb)
echo -e "Species/Population\tN_Samples\tMean_ROH_Count\tMean_Total_Length(Mb)" > $REPORT

echo "Start the concurrency controller, the current concurrency upper limit: $MAX_JOBS tasks..."

for VCF_FILE in "${VCF_FILES[@]}"; do

    POP_NAME=$(basename "$VCF_FILE" "$VCF_SUFFIX")
    echo "▶️ Submission queue: $POP_NAME ..."

    {
        # [Core Protection]: Limit each PLINK process to use up to 50,000 MB (50GB)
        MEM_MB=100000

        # Step 1: Convert VCF to PLINK high-efficiency binary format (extremely fast, memory-saving)
        $PLINK --vcf "$VCF_FILE" \
               --double-id --allow-extra-chr --chr-set 34 \
               --memory $MEM_MB \
               --make-bed \
               --out "tmp_${POP_NAME}" > /dev/null 2>&1

        # Step 2: Read the binary file and run the ROH algorithm (completed instantly)
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

        # Step Three: Check and Extract Results
        if [ -f "roh_${POP_NAME}.hom.indiv" ]; then
            # Column 4 (NSEG) is the number of ROH fragments, column 5 (KB) is the total ROH length in kb
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

            echo -e "${POP_NAME}\t${SUMMARY}" > result_line_${POP_NAME}.txt
        else
            echo -e "${POP_NAME}\t0\t0.00\t0.00" > result_line_${POP_NAME}.txt
        fi

        # Step 4: Clean the battlefield and forcibly delete the huge binary files generated in the middle to save the hard disk
        rm -f "tmp_${POP_NAME}.bed" "tmp_${POP_NAME}.bim" "tmp_${POP_NAME}.fam" "tmp_${POP_NAME}.log" "tmp_${POP_NAME}.nosex"

        echo "✅ Completed: $POP_NAME"
    } &

    # Legacy Compatible Queue Controller
    while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 5
    done

done

echo "All tasks have been issued, waiting for the last batch of calculations to be completed..."
wait

echo "All species counted! Merging reports..."
cat result_line_*.txt >> $REPORT
#rm -f result_line_*.txt

echo -e "\n========== 🎉 The analysis is done! =========="
echo "--------------------------------------------------------"
cat $REPORT | column -t
echo "--------------------------------------------------------"
echo "Detailed statistics have been saved to: $REPORT"
echo "Specific ROH fragment coordinates are saved in individual roh_species_name.hom files."
