#structure.by.admixture.sh
#!/usr/bin/bash

### input parameters
## BED file
BED_file="/DATA/03.4dTV/01.maf/01.select_snp/cycas.maffiltered.4dtv.snp.bed" ##
## input name
output_name="cycas.maffiltered.4dtv.snp" # #(Prefix of output file, required)
## sample number
sample_number=403

### software path
admixture='/TOOLS/USER/software/dongshanshan_software/admixture_linux-1.3.0/admixture';
clumpp='clumpp.pl';

## generate directory
for repeat in {11..20};do
        mkdir admixture${repeat};
done

## population structure analysis by ADMIXTURE
for K in {1..35};do
        for repeat in {11..20};do
                echo "$admixture --cv $BED_file $K -j4 -s $repeat | tee log.admixture.k${K}.${repeat}.out" > admixture.k${K}.${repeat}.sh;
        done;
done

## geneate CLUMPP input
for repeat in {11..20};do
        mv admixture.k*.${repeat}.sh admixture${repeat};
done

for repeat in {11..20};do
        for K in {1..35};do
                echo "awk '{printf NR\"\\t\"NR\"\\t\"\"(0)\"\"\\t\"4\"\\t\"\":\"\"\\t\";for(i=1;i<=NF;i++){printf \$i\"\\t\"};{printf \"\\n\"}}' admixture${repeat}/${output_name}.${K}.Q >> admixture_k${K}.indfile";
        done;
done > admixture_combine.sh

## align structure result by CLUMPP (10 is the total number of repeats, modify according to the situation)
for K in {1..35};do
        echo "perl $clumpp ${K} ${sample_number} 10 admixture";
done > clumpp.sh

#After the previous delivery task is completed, run admixture_combine.sh and clumpp.sh again, and the admixture_k*.sh script will be generated. Deliver these scripts again and get the aggregation result file admixture_k*.outfile. There is no ID in the file and it needs to be found in the .fam file.
