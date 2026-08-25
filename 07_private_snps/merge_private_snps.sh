cut -f3 sample.list.within | sort -u > poplist.txt
cat poplist.txt                              # Make sure there are no spaces in the group name
awk -F '\t' '{print $1"\t"$2 > "keep."$3".txt"}' sample.list.within

mkdir -p freq_out
for pop in $(cat poplist.txt); do
  cat > freq.${pop}.sh <<EOF
#!/bin/bash
/TOOLS/plink2 --pfile cycas.filtered \
  --keep keep.${pop}.txt --freq counts --threads 4 --out freq_out/freq_${pop}
EOF
  chmod +x freq.${pop}.sh
done
for s in freq.*.sh; do qsub -cwd -pe smp 4 "$s"; done
# Each task takes 1-3 minutes. After completion, there are 26 freq_<pop>.acount under freq_out/
