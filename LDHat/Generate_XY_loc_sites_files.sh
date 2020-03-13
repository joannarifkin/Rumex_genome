for cyt in XY
do
for LG in 7 10 3 5 8
do

vcftools --vcf /ohta/joanna.rifkin/HiCSNPData/ConvertedJosh_Felix_Pop_GQ20_Q20_Max0.5Missing_noInvariant_12-2-2019_sorted.vcf --ldhat-geno --chr L.${LG} --out RNA.${cyt}.${LG} --remove-indels --max-missing 0.8 --min-meanDP 10 --minQ 10 --max-maf 0.95 --maf 0.05 --minGQ 50 --keep ${cyt}.inds #unphased

done
done
