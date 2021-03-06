for i in "H9.NC74M" "H8.NC66F" "H10.NC84M" "H1.NC8F" "G7.NC57F" "G5.NC41M" "G4.NC33M" "G3.NC23F" "G12.NC99M" "F9.NC72M" "F8.NC64M" "F1.NC6M" "E5.NC39F" "E4.NC30F" "E3.NC21F" "E12.NC97M" "D9.NC70F" "D8.NC62M" "D1.NC4M" "C5.NC37F" "C4.NC28F" "C12.NC95M" "B9.NC68M" "B8.NC60F" "B1.NC2M" "A5.NC35F" "A4.NC26M" "A12.NC93M" "H7.NC58F" "H6.NC50M" "H5.NC42F" "H4.NC34M" "H3.NC24F" "H2.NC16F" "H12.NC100F" "H11.NC92F" "G9.NC73M" "G8.NC65M" "G6.NC49F" "G2.NC15F" "G11.NC91M" "G10.NC83M" "G1.NC7M" "F7.NC56F" "F6.NC48M" "F5.NC40F" "F4.NC31F" "F3.NC22F" "F2.NC14F" "F12.NC98F" "F11.NC90F" "F10.NC81F" "E9.NC71F" "E8.NC63F" "E7.NC55F" "E6.NC47M" "E2.NC13M" "E11.NC89M" "E10.NC80M" "E1.NC5F" "D7.NC54M" "D6.NC46M" "D5.NC38F" "D4.NC29F" "D3.NC20M" "D2.NC12M" "D12.NC96F" "D11.NC88F" "D10.NC79F" "C9.NC69F" "C8.NC61F" "C7.NC53F" "C6.NC45F" "C3.NC19F" "C2.NC11F" "C11.NC87F" "C10.NC78F" "C1.NC3M" "B7.NC52M" "B6.NC44F" "B5.NC36M" "B4.NC27M" "B3.NC18F" "B2.NC10F" "B12.NC94M" "B11.NC86F" "B10.NC77F" "A9.NC67F" "A8.NC59M" "A7.NC51M" "A6.NC43F" "A3.NC17F" "A2.NC9F" "A11.NC85F" "A10.NC75M" "A1.NC1F"

do

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  SORT_ORDER=coordinate \
      I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Pass2/11-27_Pass2_$i\Aligned.out.sam \
      O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Picard/tmp1_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	  VALIDATION_STRINGENCY=LENIENT \
     I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Picard/tmp2_$i.bam \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i


java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar MarkDuplicates \
 	  I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Picard/tmp2_$i.bam \
	  O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/NC/Picard/AnalysisReady_$i.bam \
      M=$i\marked_dup_metrics.txt \
	  REMOVE_DUPLICATES=true \
	  VALIDATION_STRINGENCY=LENIENT
	  

	  
	  

  
done	  