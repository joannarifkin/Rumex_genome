/ohta/joanna.rifkin/bcftools/bcftools mpileup --threads 20 -a AD,DP -d 3000 -Ou -f /ohta/joanna.rifkin/Genomes/Hi-C/hastate_28Sep2018_nbKXS.fasta /ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/Picard/AnalysisReady_*.bam | /ohta/joanna.rifkin/bcftools/bcftools call -vmO z -f gq -o /ohta/joanna.rifkin/HiCSNPData/Pop_RAD_seq/raw_snps_w_depth_11-12-19.vcf.gz

