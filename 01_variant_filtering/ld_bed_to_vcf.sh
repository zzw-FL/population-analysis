# Assume your file prefix is ​​mydata (i.e. mydata.bed, mydata.bim, mydata.fam)
#plink --bfile mydata --recode vcf --out output_name

# If your data is large, it is recommended to output as compressed VCF (vcf.gz)
/TOOLS/USER/software/plink/plink --bfile SAMPLE.10chr.maf.ld.snp.prune --recode vcf-fid --out SAMPLE.10chr.maf.ld.snp.prune --chr-set 34
# Or use plink2, which is faster
#plink2 --bfile mydata --recode vcf bgz --out output_name
