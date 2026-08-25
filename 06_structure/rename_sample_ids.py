#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =====================================================================
# Replace the individual number in the CLUMPP outfile with the real sample ID (directly overwrite the original file, the file name remains unchanged)
# usage:
#   python rename_ids.py ID_file [outfile ...]
#   python rename_ids.py ID_file --backup # Back up to .bak before overwriting
# ID_file format: two columns "serial number<TAB>sample ID" (the serial number starts from 1, corresponding to the outfile individual number;
#             Header rows whose first column is non-numeric will be automatically skipped; sample ID supports leading zeros such as 028)
# When outfile is not specified, all admixture_k*.outfiles in the current directory will be automatically processed.
# Note: Re-run files that have been replaced will not be processed twice (idempotent) and can be re-run with confidence.
# =====================================================================
import glob
import os
import re
import sys

NUM_RE = re.compile(r"^(\s*\d+\s+)(\d+)")   # Match the individual number in "serial number individual number" at the beginning of the line
ID_RE  = re.compile(r"^\s*(\d+)\s+(\S+)")   # ID file: "serial number sample ID"

def parse_args(argv):
    id_file = None
    backup = False
    outfiles = []
    for a in argv:
        if a in ("-b", "--backup"):
            backup = True
        elif id_file is None:
            id_file = a
        else:
            outfiles.append(a)
    return id_file, backup, outfiles

def load_ids(id_file):
    """Return (map, off): map[serial number]=sample ID; off=0 means 1-based, off=1 means 0-based"""
    mapping = {}
    maxnum = 0
    n = 0
    with open(id_file, encoding="utf-8-sig") as fh:
        for line in fh:
            m = ID_RE.match(line.strip())
            if not m:
                continue                    # Skip header/blank row
            num, name = int(m.group(1)), m.group(2)
            mapping[num] = name
            n += 1
            maxnum = max(maxnum, num)
    off = 1 if (n > 0 and maxnum == n - 1) else 0   # When 0-based, the individual number needs to be +1
    return mapping, off

def replace_ids(line, mapping, off):
    m = NUM_RE.match(line)
    if not m:
        return line                         # Non-standard rows (such as pure Q matrices), left as is
    num = int(m.group(2)) + off
    name = mapping.get(num)
    if name is None:
        return line                         # The serial number is out of range, keep the original value
    return m.group(1) + name + line[m.end():]

def main():
    id_file, backup, outfiles = parse_args(sys.argv[1:])
    if not id_file or not os.path.isfile(id_file):
        sys.exit("Usage: python rename_ids.py ID_file [--backup] [outfile ...]\n"
                 "(ID_file has two columns: serial number <TAB> sample ID; directly overwrites outfile)")
    mapping, off = load_ids(id_file)
    print(f"ID mapping: {len(mapping)} individuals, sequence number {'0-based (+1)' if off else '1-based'}")
    if not outfiles:
        outfiles = [f for f in sorted(glob.glob("admixture_k[0-9]*.outfile"))
                    if not f.endswith(".outfile.bak")]
        if not outfiles:
            sys.exit("There is no admixture_k*.outfile in the current directory")
    for f in outfiles:
        if not os.path.isfile(f):
            print(f"Skip (no file): {f}")
            continue
        tmp = f + ".tmp"
        changed = unchanged = 0
        with open(f, encoding="utf-8-sig") as fh, open(tmp, "w", encoding="utf-8", newline="\n") as fo:
            for line in fh:
                new = replace_ids(line.rstrip("\n"), mapping, off)
                if new != line.rstrip("\n"):
                    changed += 1
                else:
                    unchanged += 1
                fo.write(new + "\n")
        if backup and changed > 0:
            with open(f, encoding="utf-8-sig") as fh, open(f + ".bak", "w", encoding="utf-8", newline="\n") as fo:
                fo.write(fh.read())         # Back up the original file before replacement
        os.replace(tmp, f)                  # Atomic overwrite, file name remains unchanged
        print(f"Replaced: {f} (replace {changed} line, keep {unchanged} line)")

if __name__ == "__main__":
    main()
