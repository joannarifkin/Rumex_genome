for i in "ALBRE11" "ALBRU8" "FLJAS13" "GABEL9" "GAGLA21" "GASTA3" "LABEN5" "NCBAT4" "NCKIN2" "NCROS7" "OKBAC15" "OKRAT17" "SCBRA4" "SCMAR24" "SCPRO27" "TXATH5" "TXLIV14" "TXMTP16" "TXOAK6" "TXROS24"


do


java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  I=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Pass2/11-5_Pass2_$i\Aligned.out.sam \
	  O=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Picard/tmp1_$i.bam \
	  SORT_ORDER=coordinate \
	  VALIDATION_STRINGENCY=LENIENT



java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar MarkDuplicates \
      I=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Picard/tmp2_$i.bam \
      M=$i\marked_dup_metrics.txt \
	  REMOVE_DUPLICATES=true \
	  VALIDATION_STRINGENCY=LENIENT
	  
java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
      I=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Picard/tmp2_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/FelixHastatulusPop/AlignedtoHiC/Picard/AnalysisReady_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i

  
done	  