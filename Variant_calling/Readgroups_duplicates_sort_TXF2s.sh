for i in "A1.TX1F" "A10.TX76F" "A11.TX84M" "A12.TX93F" "A2.TX9F" "A3.TX17F" "A4.TX25M" "A5.TX33M" "A6.TX42M" "A7.TX51F" "A8.TX59F" "A9.TX68M" "B1.TX2M" "B10.TX77M" "B11.TX86F" "B12.TX94M" "B2.TX10F" "B3.TX18M" "B4.TX26M" "B5.TX34M" "B6.TX43M" "B7.TX52F" "B8.TX61M" "B9.TX69F" "C1.TX3M" "C10.TX78M" "C11.TX87F" "C12.TX95F" "C2.TX11F" "C3.TX19M" "C4.TX27M" "C5.TX35F" "C6.TX44M" "C7.TX53F" "C8.TX62F" "C9.TX70M" "D1.TX4M" "D10.TX79F" "D11.TX88M" "D12.TX96F" "D2.TX12F" "D3.TX20F" "D4.TX28M" "D5.TX36M" "D6.TX45M" "D7.TX54M" "D8.TX63F" "D9.TX71M" "E1.TX5F" "E10.TX80M" "E11.TX89F" "E12.TX97M" "E2.TX13M" "E3.TX21M" "E4.TX29M" "E5.TX37M" "E6.TX46M" "E7.TX55F" "E8.TX64F" "E9.TX72M" "F1.TX6F" "F10.TX81F" "F11.TX90M" "F12.TX98F" "F2.TX14F" "F3.TX22F" "F4.TX30F" "F5.TX38M" "F6.TX47M" "F7.TX56M" "F8.TX65M" "F9.TX73F" "G1.TX7M" "G10.TX82F" "G11.TX91F" "G12.TX99M" "G2.TX15F" "G3.TX23F" "G4.TX31M" "G5.TX39F" "G6.TX49M" "G7.TX57F" "G8.TX66F" "G9.TX74M" "H1.TX8F" "H10.TX83M" "H11.TX92M" "H12.TX100M" "H2.TX16F" "H3.TX24M" "H4.TX32M" "H5.TX40M" "H6.TX50M" "H7.TX58M" "H8.TX67M" "H9.TX75F"

do

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar SortSam \
	  SORT_ORDER=coordinate \
      I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Pass2/11-28_Pass2_$i\Aligned.out.sam \
      O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Picard/tmp1_$i.bam \
	  VALIDATION_STRINGENCY=LENIENT

java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar AddOrReplaceReadGroups \
	  VALIDATION_STRINGENCY=LENIENT \
     I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Picard/tmp1_$i.bam \
      O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Picard/tmp2_$i.bam \
	  RGID=1$i \
      RGLB=lib1$i \
      RGPL=illumina$i \
      RGPU=unit1$i \
      RGSM=sample$i



java -jar /ohta/joanna.rifkin/picard/build/libs/picard.jar MarkDuplicates \
 	  I=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Picard/tmp2_$i.bam \
	  O=/ohta/joanna.rifkin/Alignments/STAR/ComparativeTranscriptomeMap/TX/Picard/AnalysisReady_$i.bam \
      M=$i\marked_dup_metrics.txt \
	  REMOVE_DUPLICATES=true \
	  VALIDATION_STRINGENCY=LENIENT
	  

	  
	  

  
done	  