#!/bin/bash
set -e

# ==================== 1. Configuration area ====================
# [Important] Please replace this with your "filtered" high-quality VCF file name!
FINAL_VCF="/DATA/01.vcf/04.ld/10chr/SAMPLE.10chr.maf.ld.snp.prune.vcf"   

SAMPLE_LIST="/DATA/04.fis/sample.list"           # The list file you compiled (with Chinese header)
VCFTOOLS="/TOOLS/conda_down/bin/vcftools"                 # vcftools environment variable or absolute path
# ===================================================

echo "========== Phase 1: Processing sample list =========="
if [ -f "$SAMPLE_LIST" ]; then
    # Remove the first row of Chinese headers and generate a standard group mapping file
    sed '1d' $SAMPLE_LIST > pop.map
    echo "The list processing is completed and pop.map is generated, with a total of $(wc -l < ​​pop.map) sample attribution information."
else
    echo "Error: $SAMPLE_LIST file not found!"
    exit 1
fi

echo -e "\n========== Phase 2: Strict isolation calculation by species/population FIS =========="
# Create final results report file
REPORT="Final_Species_FIS_Report.txt"
echo -e "Species/Population\tN_Samples\tMean_FIS\tStatus_Eval" > $REPORT

# Get all unique species/group names
awk '{print $2}' pop.map | sort | uniq > unique_pops.txt

# Loop through each species/group
while read POP_NAME; do
    echo "Calculating independently: $POP_NAME ..."
    
    # Extract all sample names of this species
    awk -v p="$POP_NAME" '$2==p {print $1}' pop.map > temp_${POP_NAME}.txt
    N_SAM=$(wc -l < temp_${POP_NAME}.txt)
    
    # Use VCFtools to calculate independently for these specific samples (to avoid the Wahlund effect)
    $VCFTOOLS --vcf $FINAL_VCF \
              --keep temp_${POP_NAME}.txt \
              --het \
              --out res_${POP_NAME} > /dev/null 2>&1
              
    # Parse the generated .het file and calculate the average F-value of all individuals of the species
    MEAN_FIS=$(awk 'NR>1 {sum+=$5; count++} END {if(count>0) printf "%.4f", sum/count; else print "NA"}' res_${POP_NAME}.het)
    
    # Simple assessment of biological status
    if [ "$MEAN_FIS" == "NA" ]; then
        STATUS="No Data"
    else
        STATUS=$(echo "$MEAN_FIS" | awk '{
            if ($1 > 0.15) print "⚠️ High Inbreeding";
            else if ($1 > 0.05) print "⚠️ Mild Inbreeding";
            else if ($1 < -0.05) print "🌿 Potential Hybrid";
            else print "✅ Healthy (Random Mating)";
        }')
    fi
    
    # Append to the final report table
    echo -e "${POP_NAME}\t${N_SAM}\t${MEAN_FIS}\t${STATUS}" >> $REPORT
    
    # Clean temporary output files for this species
    rm temp_${POP_NAME}.txt res_${POP_NAME}.het res_${POP_NAME}.log
    
done < unique_pops.txt

echo -e "\n========== 🎉 The analysis is done! =========="
echo "The true isolation by species calculation is as follows:"
echo "--------------------------------------------------------"
cat $REPORT | column -t
echo "--------------------------------------------------------"
echo "Detailed results saved to: $REPORT"

# Clear auxiliary files
rm unique_pops.txt pop.map
