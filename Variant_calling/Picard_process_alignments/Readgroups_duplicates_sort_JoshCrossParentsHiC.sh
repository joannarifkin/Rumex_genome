for i in "TXMOM_" "TXDAD_"

do

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  SORT_ORDER=coordinate \
      I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Pass2/11-5_Pass2_$i\Aligned.out.sam \
      O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/tmp1_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	  VALIDATION_STRINGENCY=LENIENT \
     I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/tmp2_$i.bam \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i





java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar MarkDuplicates \
 	  I=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/tmp2_$i.bam \
	  O=/ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/AnalysisReady_$i.bam \
      M=$i\marked_dup_metrics.txt \
	  REMOVE_DUPLICATES=true \
	  VALIDATION_STRINGENCY=LENIENT
	  

	  
	  

  
done	  