for i in "TXMOM_" "TXDAD_"
do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R1.fastq.gz \
/ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_$i \
--runThreadN 8

done