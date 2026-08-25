#!/usr/bin/env python3
"""
get_4D_sites_ultrastrict.py — 4D  ()

():
1.  CDS  Phase  0。 Partial/EVM，！
2.  CDS Phase  GFF3 。！
3.  3 。
4. 。

:
  python3 get_4D_sites_ultrastrict.py <genome.fa> <annotation.gff3> <output.bed> <skipped.log>
"""

import sys
from collections import defaultdict

FOURFOLD_CODONS = {
    "GCT","GCC","GCA","GCG", "GGT","GGC","GGA","GGG",
    "CCT","CCC","CCA","CCG", "ACT","ACC","ACA","ACG",
    "GTT","GTC","GTA","GTG", "CTT","CTC","CTA","CTG",
    "CGT","CGC","CGA","CGG", "TCT","TCC","TCA","TCG",
}

COMPLEMENT = str.maketrans("ATGCNatgcn", "TACGNtacgn")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def parse_gff_attributes(attr_str: str) -> dict:
    attrs = {}
    for pair in attr_str.split(";"):
        pair = pair.strip()
        if "=" in pair:
            k, v = pair.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs

def load_genome(fasta_path: str) -> dict:
    genome = {}
    current_chrom = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom is not None:
                    genome[current_chrom] = "".join(current_seq)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_chrom is not None:
        genome[current_chrom] = "".join(current_seq)
    return genome

def load_cds_from_gff(gff_path: str) -> list:
    cdss = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "CDS":
                continue
            attrs = parse_gff_attributes(cols[8])
            parent = attrs.get("Parent", "")
            if not parent:
                continue
            phase_str = cols[7]
            phase = int(phase_str) if phase_str in ("0", "1", "2") else 0
            cdss.append({
                "chrom":  cols[0],
                "start":  int(cols[3]),
                "end":    int(cols[4]),
                "strand": cols[6],
                "phase":  phase,
                "parent": parent,
            })
    return cdss

def process_transcript(genome, cdss, log_fh):
    if not cdss:
        return set()

    strand = cdss[0]["strand"]
    chrom = cdss[0]["chrom"]
    parent_id = cdss[0]["parent"]

    if chrom not in genome:
        return set()

    if strand == "+":
        cdss_sorted = sorted(cdss, key=lambda c: c["start"])
    else:
        cdss_sorted = sorted(cdss, key=lambda c: c["start"], reverse=True)

    # ========================================================
    # QC 1: The first CDS phase must be 0 (no Partial Gene is tolerated)
    # ========================================================
    if cdss_sorted[0]["phase"] != 0:
        log_fh.write(f"[SKIP_PHASE_0] {parent_id}: first CDS phase={cdss_sorted[0]['phase']} (determined as incomplete or EVM modified model, discarded)\n")
        return set()

    # ========================================================
    # Quality Control 2: Internal Phases must be perfectly connected (experts’ ultimate advice)
    # ========================================================
    for i in range(len(cdss_sorted) - 1):
        curr = cdss_sorted[i]
        nxt = cdss_sorted[i + 1]
        curr_len = curr["end"] - curr["start"] + 1
        
        # curr is the number of missing codon bases left for nxt
        leftover = (curr_len - curr["phase"]) % 3
        # nxt How many bases need to be skipped to open a new codon
        expected_next_phase = (3 - leftover) % 3
        
        if nxt["phase"] != expected_next_phase:
            log_fh.write(f"[SKIP_INTERNAL_PHASE] {parent_id}: Exon {i+1} and Exon {i+2} Phase connection exception (expected {expected_next_phase}, actual {nxt['phase']}). There is a frameshift internally, which is discarded. \n")
            return set()

    # ---------- Splicing sequence and mapping coordinates ----------
    cds_seq = ""
    genome_pos = []

    if strand == "+":
        for cds in cdss_sorted:
            gseq = genome[chrom][cds["start"] - 1 : cds["end"]]
            if not gseq: continue
            cds_seq += gseq
            for offset_in_cds in range(len(gseq)):
                genome_pos.append(cds["start"] + offset_in_cds)
    else:
        for cds in cdss_sorted:
            gseq = genome[chrom][cds["start"] - 1 : cds["end"]]
            if not gseq: continue
            rc = revcomp(gseq)
            cds_seq += rc
            for offset_in_rc in range(len(rc)):
                genome_pos.append(cds["end"] - offset_in_rc)

    # ========================================================
    # Quality control 3: The total splicing length must be a multiple of 3
    # ========================================================
    if len(cds_seq) % 3 != 0:
        log_fh.write(f"[SKIP_LENGTH] {parent_id}: CDS is perfectly connected but the total length ({len(cds_seq)} bp) is not divisible by 3 (suspected to be end truncation). throw away. \n")
        return set()

    if len(cds_seq) < 3:
        return set()

    # ---------- Scan 4D sites codon by codon ----------
    sites = set()
    for i in range(0, len(cds_seq), 3):
        codon = cds_seq[i:i + 3]
        if codon in FOURFOLD_CODONS:
            pos3_1based = genome_pos[i + 2]
            sites.add((chrom, pos3_1based - 1, pos3_1based))

    return sites


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 get_4D_sites_ultrastrict.py <genome.fa> <annotation.gff3> <output.bed> <skipped.log>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    gff_path = sys.argv[2]
    out_path = sys.argv[3]
    log_path = sys.argv[4]

    print("Loading reference genome...")
    genome = load_genome(fasta_path)
    
    print("Parsing GFF3 CDS comments...")
    all_cdss = load_cds_from_gff(gff_path)

    transcripts = defaultdict(list)
    for cds in all_cdss:
        transcripts[cds["parent"]].append(cds)
    
    print(f"Start scanning 4D sites of {len(transcripts)} transcripts...")
    print(f"(Perform ultimate quality control: first phase verification + internal phase connection verification + multiples of 3 verification)")
          
    all_sites = set()
    tx_with_4d = set()

    with open(log_path, "w") as log_fh:
        log_fh.write("========== Log of transcripts discarded by Ultimate Quality Control ==========\n")
        
        for parent, cdss in transcripts.items():
            sites = process_transcript(genome, cdss, log_fh)
            if sites:
                tx_with_4d.add(parent)
                all_sites.update(sites)

    # Write deduplication BED
    sorted_sites = sorted(all_sites, key=lambda x: (x[0], x[1]))
    with open(out_path, "w") as out_fh:
        for chrom, start, end in sorted_sites:
            out_fh.write(f"{chrom}\t{start}\t{end}\n")

    print(f"\n✅ Perfect ending! Among {len(tx_with_4d)} absolutely perfect transcripts, {len(all_sites)} 4D sites were extracted!")
    print(f"➜ Output BED: {out_path}")
    print(f"➜ Quality control log: {log_path}")

if __name__ == "__main__":
    main()
