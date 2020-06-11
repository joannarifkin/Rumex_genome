./BuildDatabase -name HastMale -engine ncbi filtered_5000_a.lines.fasta 
nohup ./RepeatModeler -engine ncbi -pa 3 -database  Dovetail -recoverDir RM_4574.ThuDec211156042017 >& recoverrun1221Dovetail.out

#Concatenate all repeat libraries
cat BuckwheatFull-families.fa Dovetail-families.fa HastMale-families.fa HastFemale-families.fa > RumexAndBuckwheatRepeats.fa
#Make classified library from repeat libraries
RepeatClassifier -consensi RumexAndBuckwheatRepeats.fa 

####filter assembly####
	##remove short contigs | #remove complex symbols in contig names (; & =) & all to uppercase
	#/ohta/apps/kentUtils/bin/faSize -detailed /ohta/felix.beaudry/assemblies/hastGenome/dovetail/hastate_28Sep2018_nbKXS.fasta
	perl removesmalls.pl 50000 hastate_28Sep2018_nbKXS.fasta | sed 's/=/-/' | sed 's/;/-/' | awk '{print toupper($0)}' > hastate_28Sep2018_nbKXS_50KUP.fasta

	##softmask repeats (RepeatMaskerCommands.txt) based on buckwheat + rumex library
	RepeatMasker -pa 4 -xsmall -dir soft -lib RumexAndBuckwheatRepeats.fa.classified hastate_28Sep2018_nbKXS_50KUP.fasta  

	##make alignment database
	mkdir hastate_28Sep2018_nbKXS_50KUP
	STAR --runMode genomeGenerate --genomeDir hastate_28Sep2018_nbKXS_50KUP --genomeFastaFiles hastate_28Sep2018_nbKXS_50KUP.masked.fasta --runThreadN 10
	java -jar picard.jar CreateSequenceDictionary R=hastate_28Sep2018_nbKXS_50KUP.masked.fasta O=hastate_28Sep2018_nbKXS_50KUP.masked.dict

	####BRAKER lignments####
		##female leaf
		STAR --genomeDir hastate_28Sep2018_nbKXS_50KUP --readFilesIn 1.TF1_R1_clean.fastq.gz 1.TF1_R2_clean.fastq.gz --twopassMode Basic --outFileNamePrefix leaf_TX1F. --runThreadN 10 --readFilesCommand zcat

		##male leaf
		STAR --genomeDir hastate_28Sep2018_nbKXS_50KUP --readFilesIn TXROS24_R1.fastq TXROS24_R2.fastq --twopassMode Basic --outFileNamePrefix leaf_TXROS24M. --runThreadN 10

		for ind in 1 2
		do
		STAR --genomeDir hastate_28Sep2018_nbKXS_50KUP --readFilesIn Hast_TX${ind}B_R1_clean.fastq  Hast_TX${ind}B_R2_clean.fastq --twopassMode Basic --outFileNamePrefix pollen_TX${ind}B. --runThreadN 10 #--readFilesCommand zcat
		done

		for ind in 17_TXMale 24_TXFem
		do
		TAR --genomeDir hastate_28Sep2018_nbKXS_50KUP --readFilesIn ${ind}Flower_R1_clean.fastq.gz ${ind}Flower_R1_clean.fastq.gz --twopassMode Basic --outFileNamePrefix flower_${ind}. --runThreadN 10 --readFilesCommand zcat
		done

		for ind in  pollen_TX1B flower_17_TXMale flower_24_TXFem pollen_TX2B leaf_TXROS24M leaf_TX1F
		do
		echo "cleaning ${ind}"
		##input BAM file must be sorted by target (=genome) sequence names and within the sequences by begin coordinates
		java -jar picard.jar SortSam SORT_ORDER=coordinate I=${ind}.Aligned.out.sam O=${ind}.sort.temp.bam VALIDATION_STRINGENCY=LENIENT

		java -jar picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT I=${ind}.sort.temp.bam O=${ind}.add.bam RGID=${ind} RGLB=RNAseq RGPL=RNAseq RGPU=unit1 RGSM=${ind}

		java -jar picard.jar MarkDuplicates I=${ind}.add.temp.bam O=${ind}.ddpl.bam M=${ind}.ddpl_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

		done

	####BRAKER####
	braker.pl --genome=hastate_28Sep2018_nbKXS_50KUP.masked.fasta --species=Rhast --softmasking --gff3 --cores=8 --bam=pollen_TX1B.add.bam,flower_17_TXMale.add.bam,flower_24_TXFem.add.bam,pollen_TX2B.add.bam --PYTHON3_PATH=/usr/bin --AUGUSTUS_CONFIG_PATH=./Augustus/config/ --augustus_args=--species=arabidopsis --rounds=4


	####gene filter####
	##make assembly that is just genes
		awk '$3 ~ "gene" {print}' augustus.hints.gff3 | awk '$2 ~ "AUGUSTUS" {print}' | awk '$6 < 1 {print}'  |  sed 's/=//' | sed 's/;//' | awk '{print $1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >augustus.genes.bed
		bedtools getfasta -fi hastate_28Sep2018_nbKXS_50KUP.masked.fasta -bed augustus.genes.bed -fo transcriptomePrediction.fa -s -name 
		#fold -w 60 transcriptomePrediction.long.fa >transcriptomePrediction.fa
		awk '$1 ~ ">" {print}' transcriptomePrediction.fa | awk 'split($1,a,">") {print a[2]}' >base.txt

		mkdir transcriptomePrediction
		STAR --runMode genomeGenerate --genomeDir transcriptomePrediction --genomeFastaFiles transcriptomePrediction.fa --runThreadN 10 --limitGenomeGenerateRAM 60000000000

		##alignment for coverage analysis##
		##RNA##
		STAR --genomeDir transcriptomePrediction --readFilesIn 1.TF1_R1_clean.fastq.gz pop/1.TF1_R2_clean.fastq.gz --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix leaf_TX1F. --runThreadN 10 --readFilesCommand zcat
		STAR --genomeDir transcriptomePrediction --readFilesIn TXROS24_R1.fastq TXROS24_R2.fastq --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix leaf_TXROS24M. --runThreadN 10

		for ind in 1 2
		do
		STAR --genomeDir transcriptomePrediction --readFilesIn Hast_TX${ind}B_R1_clean.fastq  Hast_TX${ind}B_R2_clean.fastq --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix pollen_TX${ind}B. --runThreadN 10 #--readFilesCommand zcat
		done

		for ind in 17_TXMale 24_TXFem
		do
		STAR --genomeDir transcriptomePrediction --readFilesIn ${ind}Flower_R1_clean.fastq.gz ${ind}Flower_R2_clean.fastq.gz --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix flower_${ind}. --runThreadN 10 --readFilesCommand zcat 
		done

		for ind in  flower_17_TXMale flower_24_TXFem #pollen_TX1B pollen_TX2B leaf_TXROS24M leaf_TX1F
		do
		samtools index ${ind}.Aligned.sortedByCoord.out.bam
		samtools idxstats ${ind}.Aligned.sortedByCoord.out.bam >${ind}.rna.reads.txt
		done

		##genomic##
		rm transcriptomePrediction.fa-enc.2.ngm transcriptomePrediction.fa-ht-13-2.3.ngm
		ngm -t 10 -r transcriptomePrediction.fa -b -p -1 TxOAK20m_R1_clean.fastq.gz -2 TxOAK20m_R2_clean.fastq.gz -o TXOAK20m.Aligned.out.bam 
		ngm -t 10 -r transcriptomePrediction.fa -b -p -1 HI.3472.002.D709---D502.r_hastatulus_Tx_PE_female_R1_clean.fastq.gz -2 HI.3472.002.D709---D502.r_hastatulus_Tx_PE_female_R2_clean.fastq.gz -o TXROS2f.Aligned.out.bam
		ngm -t 10 -r transcriptomePrediction.fa -b -p -1 HI.2651.002.Index_2.Rumex_female_R1_clean.fastq.gz -2 HI.2651.002.Index_2.Rumex_female_R2_clean.fastq.gz -o SCMAR2f.Aligned.out.bam
		ngm -t 10 -r transcriptomePrediction.fa -b -p -1 HI.2651.001.Index_13.Rumex_male_R1_clean.fastq.gz -2 HI.2651.001.Index_13.Rumex_male_R2_clean.fastq.gz -o SCMAR1m.Aligned.out.bam
		java -jar picard.jar CreateSequenceDictionary R=transcriptomePrediction.fa O=transcriptomePrediction.dict


		for ind in TXOAK20m TXROS2f SCMAR2f SCMAR1m
		do
		java -jar picard.jar SortSam SORT_ORDER=coordinate I=${ind}.Aligned.out.bam O=${ind}.Aligned.sortedByCoord.out.bam VALIDATION_STRINGENCY=LENIENT
		samtools index ${ind}.Aligned.sortedByCoord.out.bam
		samtools idxstats ${ind}.Aligned.sortedByCoord.out.bam >${ind}.gen.reads.txt
		done

	####annotation####
	#goodgenes filtered in R, annotationFilter.R #
	awk 'split($1,a,"ID") {print a[2]}' goodgenes.list >goodgenes.ls

	##filter##
	rm goodgenes.gtf
	while read loc
	do
	echo "looking for ${loc}"
	awk -v var=${loc} '$3 == "gene" && $9 == var {print}' augustus.hints.gtf >>goodgenes.gtf
	awk -v var=${loc}.t '$3 == "transcript" && $9 ~ var {print}' augustus.hints.gtf >>goodgenes.gtf
	awk -v var=${loc}.t ' $10 ~ var {print}' augustus.hints.gtf >>goodgenes.gtf
	done <goodgenes.ls

	sed 's/SCNBKXS/ScnbKXS/' goodgenes.gff | sed 's/-HRSCAF-/_HRSCAF_/' >goodGenes.gff
	python3 gff_Interactive_Chromonomer_LG_position_switcher_FB.py goodGenes.gff ann2/goodGenes.conv.chrom.gff
	awk '{print $4"\t"$7"\t"$8"\t"$5"\t"$6"\t"$9"\t"$10"\t"$11"\t"$12}' goodGenes.conv.chrom.gff >goodGenes.chrom.gff 
	rm goodGenes.conv.chrom.gff goodGenes.gff

	awk '$3 ~ "gene"  {print}' goodGenes.chrom.gff | awk ' split($9,a,"=") split(a[2],b,";")  {print "ID"b[1]"\t"$1"\t"$4"\t"$5}' >goodGenes.chrom.pos

	awk '$3 ~ "gene" {print}' goodgenes.gff | awk '$2 ~ "AUGUSTUS" {print}' | awk '$6 < 1 {print}'  |  sed 's/=//' | sed 's/;//' | awk '{print $1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' >goodgenes.bed
	./bedtools2/bin/fastaFromBed -fi hastate_28Sep2018_nbKXS_50KUP.masked.fasta -bed goodgenes.bed -fo goodGenes.fa -s -name 

	#BUSCO run locally
	docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.5_cv1 busco -i goodGenes.fa -o goodGenes.busco -m tran -l viridiplantae_odb10 -f

docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v4.0.5_cv1 busco -i goodGenes.fa -o goodGenes.busco -m tran -l eukaryota_odb10 -f
