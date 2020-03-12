for i in "8_Index-5.NCMOM_" "7_Index-11.NCDAD_"

do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R1_clean.fastq.gz \
/ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R2_clean.fastq.gz \
--readFilesCommand zcat \
--sjdbFileChrStartEnd \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_TXMOM_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_TXDAD_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_8_Index-5.NCMOM_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_7_Index-11.NCDAD_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_27.NCF1_M6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_25.NCF1_M5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_23.NCF1_M4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_22.NCF1_M3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_21.NCF1_M2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_20.NCF1_M1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_19.NCF1_F6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_18.NCF1_F5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_16.NCF1_F4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_15.NCF1_F3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_14.NCF1_F2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_2_Index_13.NCF1_F1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-9.TXF1_M3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-8.TXF1_M2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-7.TXF1_M1_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-6.TXF1_F6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-5.TXF1_F5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-4.TXF1_F4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-3.TXF1_F3_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-2.TXF1_F2_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-12.TXF1_M6_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-11.TXF1_M5_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-10.TXF1_M4_SJ.out.tab \
/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_1_Index-1.TXF1_F1_SJ.out.tab \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass2/11-5_Pass2_$i \
--runThreadN 16

done
