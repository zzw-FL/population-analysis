# Add -upDownStreamLen 2000 to strictly limit the upstream and downstream areas
/TOOLS/conda_down/bin/java -Xmx124g -jar /TOOLS/snpEff/snpEff.jar -v   -upDownStreamLen 2000 -s SAMPLE_snpEff_summary.html REFERENCE SAMPLE.10chr.maf.ld.snp.prune.vcf > SAMPLE.ann.vcf
