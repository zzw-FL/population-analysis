/TOOLS/USER/software/sourcesoftware/bcftools/bcftools view -S /DATA/01.vcf/POPULATION.list -Oz -o POPULATION.list.ld.snp.vcf.gz /DATA/01.vcf/04.ld/10chr/SAMPLE.10chr.maf.ld.snp.prune.vcf
/TOOLS/USER/software/sourcesoftware/bcftools/bcftools index -t POPULATION.list.ld.snp.vcf.gz
