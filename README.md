# Population analysis

Curated population-genomics workflows from the cluster `09.china_cycas`
project. Representative scripts were selected and sample-specific or
function-duplicated copies were removed.

## Directories

| Directory | Function |
| --- | --- |
| `01_variant_filtering/` | MAF/LD filtering, sample subsetting and ID renaming |
| `02_phylogenetic_tree/` | VCF to PHYLIP conversion, RAxML-NG and FastTree |
| `03_fourfold_sites/` | 4DTv site extraction and strict filtering |
| `04_fis/` | Population inbreeding coefficient estimation |
| `05_roh/` | Runs of homozygosity calculation and summarization |
| `06_structure/` | ADMIXTURE, CLUMPP, CV plotting and structure plots |
| `07_private_snps/` | Private-SNP frequency, filtering and summarization |

## Privacy notes

- Absolute cluster paths, user names, project identifiers and internal sample
  names were replaced with placeholders such as `/DATA`, `/TOOLS`, `/TOOLS/USER`,
  `USER`, `PROJECT`, `SAMPLE`, and `POPULATION`.
- Raw VCF/FASTQ/BAM data, sample metadata, result tables and figures were not
  uploaded.

## License

See [LICENSE](LICENSE).
