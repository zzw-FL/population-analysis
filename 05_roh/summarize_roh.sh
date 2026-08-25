#!/bin/bash
set -e

echo "========== Start parsing PLINK ROH results =========="

# Prepare three output tables
IND_OUT="Table1_ROH_Individual_Summary.txt"
POP_OUT="Table2_ROH_Population_Summary.txt"
SEG_OUT="Table3_ROH_Segment_Distribution.txt"

# Write the table header (for the convenience of subsequent drawing and reading, the length is uniformly converted into Mb)
echo -e "Population\tSample\tN_ROH\tAverage_Length_Mb\tTotal_Length_Mb" > $IND_OUT
echo -e "Population\tN_Samples\tMean_N_ROH\tMean_Total_Length_Mb\tPopulation_Avg_Segment_Mb" > $POP_OUT
echo -e "Population\tSample\tROH_Length_Mb" > $SEG_OUT

# Find all .hom files in the current directory
shopt -s nullglob
HOM_FILES=(roh_*.hom)
if [ ${#HOM_FILES[@]} -eq 0 ]; then
    echo "❌ Error: roh_*.hom files not found in current directory! Please confirm whether the previous script ran successfully."
    exit 1
fi

echo "✅ Found ROH data for ${#HOM_FILES[@]} species; analyzing and summarizing across multiple dimensions..."

for HOM_FILE in "${HOM_FILES[@]}"; do
    # Extract group/species name
    POP_NAME=$(basename "$HOM_FILE" .hom | sed 's/^roh_//')
    
    INDIV_FILE="roh_${POP_NAME}.hom.indiv"
    
    if [ ! -f "$INDIV_FILE" ]; then
        echo "⚠️ WARNING: $INDIV_FILE not found, $POP_NAME skipped..."
        continue
    fi

    # ---------------------------------------------------------
    # 1. Extract segment length distribution (Segment-level)
    # Column 2 of the .hom file is the sample name (IID) and column 9 is the length (KB)
    # ---------------------------------------------------------
    awk -v pop="$POP_NAME" 'NR>1 {
        print pop "\t" $2 "\t" $9/1000
    }' "$HOM_FILE" >> $SEG_OUT

    # ---------------------------------------------------------
    # 2. Individual-level statistics
    # .hom.indiv file sample name in column 2, ROH quantity (NSEG) in column 4, total length (KB) in column 5
    # ---------------------------------------------------------
    awk -v pop="$POP_NAME" 'NR>1 {
        nseg = $4
        total_kb = $5
        # Prevent error reporting when dividing by 0
        avg_mb = (nseg > 0) ? (total_kb/1000)/nseg : 0
        print pop "\t" $2 "\t" nseg "\t" avg_mb "\t" total_kb/1000
    }' "$INDIV_FILE" >> $IND_OUT

    # ---------------------------------------------------------
    # 3. Population-level summary statistics
    # ---------------------------------------------------------
    awk -v pop="$POP_NAME" 'NR>1 {
        nseg += $4
        total_kb += $5
        count++
    } END {
        if(count > 0) {
            mean_nseg = nseg / count
            mean_total_mb = (total_kb / 1000) / count
            # Average fragment length of the population = total length of all fragments / total number of all fragments
            pop_avg_seg = (nseg > 0) ? (total_kb / 1000) / nseg : 0
            
            printf "%s\t%d\t%.2f\t%.2f\t%.2f\n", pop, count, mean_nseg, mean_total_mb, pop_avg_seg
        }
    }' "$INDIV_FILE" >> $POP_OUT

done

echo -e "\n========== 🎉 Analysis completed! =========="
echo "Three statistical reports are generated that perfectly match your needs:"
echo "1. Individual level statistics table: $IND_OUT"
echo "2. Population level summary table: $POP_OUT"
echo "3. Fragment distribution table: $SEG_OUT (can be used to draw violin plots/histograms in R language)"
