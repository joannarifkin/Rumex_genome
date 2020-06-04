while read i
do



java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  I=/ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/RawBam/$i.bam \
	  O=/ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/Picard/tmp1_$i.bam \
	  SORT_ORDER=coordinate \
	  VALIDATION_STRINGENCY=LENIENT

	  
java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
      I=/ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/NGM/Hastatulus/HiC/PopRadSeq/Picard/AnalysisReady_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i

  

done < /ohta/joanna.rifkin/Alignments/NGM/Hastatulus/ChicagoAssembly/PopRadSeq/FelixRADSeqPopList.txt