/TOOLS/USER/software/sourcesoftware/bcftools/bcftools view -S /DATA/01.vcf/POPULATION.list -Oz -o POPULATION.list.maf.snp.vcf.gz /DATA/01.vcf/03.maf.merge/SAMPLE.10chr.maf.snp.vcf.gz
/TOOLS/USER/software/sourcesoftware/bcftools/bcftools index -t POPULATION.list.maf.snp.vcf.gz
