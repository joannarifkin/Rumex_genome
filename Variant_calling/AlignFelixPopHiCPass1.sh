for i in "ALBRE11" "ALBRU8" "FLJAS13" "GABEL9" "GAGLA21" "GASTA3" "LABEN5" "NCBAT4" "NCKIN2" "NCROS7" "OKBAC15" "OKRAT17" "SCBRA4" "SCMAR24" "SCPRO27" "TXATH5" "TXLIV14" "TXMTP16" "TXOAK6" "TXROS24" "FLHAM4" "NCELI16" "NCHIC16"
do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R1.fastq \
/ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R2.fastq \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_$i \
--runThreadN 16

done
