for i in "FLHAM4" "NCELI16" "NCHIC16"


do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R1.fastq \
/ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R2.fastq \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_$i \
--runThreadN 30

done