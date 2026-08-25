#!/bin/bash
set -e

# ==================== Configuration area (please modify here) ====================
ANN_VCF="SAMPLE.ann.vcf"               # The VCF file you just annotated with SnpEff
REF_FASTA="REFERENCE.fasta"    # Your reference genome file (name changed to the purely digital version)
SNP_SIFT_JAR="/TOOLS/snpEff/SnpSift.jar" # Absolute path to SnpSift
# =============================================================

echo "========== Phase 1: Extract synonymous mutations (Synonymous) using SnpSift =========="
/TOOLS/conda_down/bin/java -jar $SNP_SIFT_JAR filter "ANN[*].EFFECT has 'synonymous_variant'" $ANN_VCF > step1_synonymous.vcf

echo "========== Stage 2: Extract transversion sites (Tv) using bcftools =========="
# Only variations in purine and pyrimidine exchanges are retained
/TOOLS/USER/software/sourcesoftware/bcftools/bcftools view -i '((REF="A" && ALT="C") || (REF="A" && ALT="T") || (REF="C" && ALT="A") || (REF="C" && ALT="G") || (REF="G" && ALT="C") || (REF="G" && ALT="T") || (REF="T" && ALT="A") || (REF="T" && ALT="G"))' step1_synonymous.vcf > step2_syn_Tv.vcf

echo "========== Stage 3: Extract tight 4DTv using Python Exact Alignment Fasta =========="
# FASTA must be indexed first
if [ ! -f "${REF_FASTA}.fai" ]; then
    echo "Building samtools index for reference genome..."
    /TOOLS/software/bin/samtools faidx $REF_FASTA
fi

# Run Python script for final refinement
/TOOLS/software/bin/python3 extract_strict_4DTv.py step2_syn_Tv.vcf final_strict_4DTv.vcf $REF_FASTA

echo "========== Stage 4: Clean temporary files and output final results =========="
rm step1_synonymous.vcf step2_syn_Tv.vcf
echo "Congratulations! The purest 4DTv variant file has been generated: final_strict_4DTv.vcf"
