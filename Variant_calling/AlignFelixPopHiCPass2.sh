for i in "ALBRE11" "ALBRU8" "FLJAS13" "GABEL9" "GAGLA21" "GASTA3" "LABEN5" "NCBAT4" "NCKIN2" "NCROS7" "OKBAC15" "OKRAT17" "SCBRA4" "SCMAR24" "SCPRO27" "TXATH5" "TXLIV14" "TXMTP16" "TXOAK6" "TXROS24" "FLHAM4" "NCELI16" "NCHIC16"
do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R1.fastq \
/ohta/felix.beaudry/rawSequence/RNAseq/pop/$i\_R2.fastq \
--sjdbFileChrStartEnd \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_TXROS24SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_TXOAK6SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_TXMTP16SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_TXLIV14SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_TXATH5SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_SCPRO27SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_SCMAR24SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_SCBRA4SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_OKRAT17SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_OKBAC15SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_NCROS7SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_NCKIN2SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_NCBAT4SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_LABEN5SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_GASTA3SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_GAGLA21SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_GABEL9SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_FLJAS13SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_ALBRU8SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_ALBRE11SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_FLHAM4SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_NCELI16SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass1/11-2_Pass1_NCHIC16SJ.out.tab \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass2/11-5_Pass2_$i \
--runThreadN 16

done
