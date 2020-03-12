for i in "10.TM4_" "11.TM5_" "12.TM6_" "1.TF1_" "2.TF2_" "3.TF3_" "4.TF4_" "5.TF5_" "6.TF6_" "7.TM1_" "8.TM2_" "9.TM3_" "13.NF1_" "14.NF2_" "15.NF3_" "16.NF4_" "18.NF5_" "19.NF6_" "20.NM1_" "21.NM2_" "22.NM3_" "23.NM4_" "25.NM5_" "27.NM6_"


do

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Pass2/11-5_Pass2_$i\Aligned.out.sam \
	  O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Picard/tmp1_$i.bam \
	  SORT_ORDER=coordinate \
	  VALIDATION_STRINGENCY=LENIENT



java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar MarkDuplicates \
      I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Picard/tmp2_$i.bam \
      M=$i\marked_dup_metrics.txt \
	  REMOVE_DUPLICATES=true \
	  VALIDATION_STRINGENCY=LENIENT
	  
java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
      I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Picard/tmp2_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusPop/Aligned_To_HiC_Assembly/Picard/AnalysisReady_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i

  
done	  