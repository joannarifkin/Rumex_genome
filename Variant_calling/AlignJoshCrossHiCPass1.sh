for i in "8_Index-5.NCMOM_" "7_Index-11.NCDAD_" "6_Index-4.TXM_" "6_Index-2.TX3_" "6_Index-1.TX1_" "5_Index-4.TXM_" "5_Index-2.TX3_" "5_Index-1.TX1_" "2_Index_27.NCF1_M6_" "2_Index_25.NCF1_M5_" "2_Index_23.NCF1_M4_" "2_Index_22.NCF1_M3_" "2_Index_21.NCF1_M2_" "2_Index_20.NCF1_M1_" "2_Index_19.NCF1_F6_" "2_Index_18.NCF1_F5_" "2_Index_16.NCF1_F4_" "2_Index_15.NCF1_F3_" "2_Index_14.NCF1_F2_" "2_Index_13.NCF1_F1_" "1_Index-9.TXF1_M3_" "1_Index-8.TXF1_M2_" "1_Index-7.TXF1_M1_" "1_Index-6.TXF1_F6_" "1_Index-5.TXF1_F5_" "1_Index-4.TXF1_F4_" "1_Index-3.TXF1_F3_" "1_Index-2.TXF1_F2_" "1_Index-12.TXF1_M6_" "1_Index-11.TXF1_M5_" "1_Index-10.TXF1_M4_" "1_Index-1.TXF1_F1_"

do

STAR --genomeDir /ohta/joanna.rifkin/Genomes/Hi-C \
--readFilesIn /ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R1_clean.fastq.gz \
/ohta/felix.beaudry/rawSequence/RNAseq/josh/cross/$i\R2_clean.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix /ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass1/11-1_Pass1_$i \
--runThreadN 16

done
