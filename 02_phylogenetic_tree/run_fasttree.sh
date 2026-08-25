#!/bin/bash
#SBATCH -J fasttree_cycas
#SBATCH -c 48
#SBATCH --mem=32G
#SBATCH -o fasttree.log
#SBATCH -e fasttree.err

FASTTREE=/TOOLS/fasttree/FastTree

INPUT=SAMPLE.10chr.maf.ld.snp.prune.min4.phy.varsites.fa
OUT=cycas_fasttree.nwk

echo "Start FastTree..."

$FASTTREE \
    -nt \
    -gtr \
    -cat 20 \
    -log cycas_fasttree.log \
    $INPUT \
    > $OUT

echo "Done!"
