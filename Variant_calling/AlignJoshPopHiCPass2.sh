for i in "10.TM4_" "11.TM5_" "12.TM6_" "1.TF1_" "2.TF2_" "3.TF3_" "4.TF4_" "5.TF5_" "6.TF6_" "7.TM1_" "8.TM2_" "9.TM3_" "13.NF1_" "14.NF2_" "15.NF3_" "16.NF4_" "18.NF5_" "19.NF6_" "20.NM1_" "21.NM2_" "22.NM3_" "23.NM4_" "25.NM5_" "27.NM6_"

do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/josh/pop/$i\R1_clean.fastq \
/ohta/felix.beaudry/rawSequence/RNAseq/josh/pop/$i\R2_clean.fastq \
--sjdbFileChrStartEnd \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_9.TM3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_8.TM2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_7.TM1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_6.TF6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_5.TF5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_4.TF4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_3.TF3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_27.NM6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_25.NM5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_23.NM4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_22.NM3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_21.NM2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_20.NM1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2.TF2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_19.NF6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_18.NF5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_16.NF4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_15.NF3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_14.NF2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_13.NF1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_12.TM6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_11.TM5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_10.TM4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1.TF1_SJ.out.tab \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass2/11-5_Pass2_$i \
--runThreadN 16

done
