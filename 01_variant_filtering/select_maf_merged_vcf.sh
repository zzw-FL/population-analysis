/TOOLS/USER/software/sourcesoftware/bcftools/bcftools view -S /DATA/02.cycas/03.filter/pop25/POPULATION/POPULATION.txt -Oz -o POPULATION.xx.allhard.snp.vcf.gz SAMPLE.10chr.maf.snp.vcf.gz
/TOOLS/USER/software/sourcesoftware/bcftools/bcftools index -t POPULATION.xx.allhard.snp.vcf.gz
