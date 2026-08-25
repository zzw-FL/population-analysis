#!/usr/bin/env python3
import sys
from pysam import VariantFile, FastaFile

if len(sys.argv) != 4:
    print("Usage: python3 extract_strict_4DTv.py <input.vcf> <output.vcf> <reference.fasta>")
    sys.exit(1)

file_in = sys.argv[1]
file_out = sys.argv[2]
fafile = sys.argv[3]

# Strict 4-fold degenerate amino acid first two codon combinations (positive strand)
codon = set(["TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG"])
rev_dict = dict(A='T', T='A', C='G', G='C', N='N')

bcf_in = VariantFile(file_in)
bcf_out = VariantFile(file_out, "w", header=bcf_in.header)
fa_in = FastaFile(fafile)

kept_count = 0
dropped_count = 0

for rec in bcf_in.fetch():
    # Make sure there is an ANN field (note for SnpEff)
    if 'ANN' not in rec.info:
        continue
    
    # Extract SnpEff annotation information (take the first most important impact)
    info = rec.info['ANN'][0].split('|')
    
    # Filter layer 1: must be synonymous mutation
    if "synonymous_variant" not in info[1]:
        continue
        
    # Filter layer 2: Must be at position 3 of the codon (parse HGVS.c information such as c.123A>G)
    hgvs_c = info[9]
    if not hgvs_c.startswith('c.'):
        continue
        
    try:
        # Extract variant cDNA coordinate numbers
        pos_num_str = ''.join([c for c in hgvs_c[2:] if c.isdigit()])
        if int(pos_num_str) % 3 != 0:
            continue
    except ValueError:
        continue

    # Filtering layer 3: Verify whether it is a true 4D codon based on the reference genome
    try:
        # Extract the reference base before mutation (REF)
        ref_base = hgvs_c[-3] if '>' in hgvs_c else rec.ref 
        
        # If the gene where the mutation is located is on the positive strand (VCF REF base is consistent with the HGVS REF base annotated by SnpEff)
        if rec.ref == ref_base:
            # Extract the first 2 bases of the mutation site
            pre = fa_in.fetch(rec.chrom, rec.pos - 3, rec.pos - 1).upper()
        # If the gene where the mutation is located is on the negative strand
        else:
            # Extract the 2 bases behind the mutation site and take the reverse complement
            tmp = fa_in.fetch(rec.chrom, rec.pos, rec.pos + 2).upper()
            pre = rev_dict.get(tmp[1], 'N') + rev_dict.get(tmp[0], 'N')
            
        if pre in codon:
            bcf_out.write(rec)
            kept_count += 1
        else:
            dropped_count += 1
            
    except Exception as e:
        continue

print(f"Tight filtering is completed: {kept_count} real 4DTv sites are retained, and {dropped_count} non-4D pseudo sites are eliminated.")
