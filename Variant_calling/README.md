#### Variant_calling

	#### Align
		Index_genome.txt #Script to index our draft genome assembly for use with STAR
=======
 
Index_genome.txt #Script to index our draft genome assembly for use with STAR (ngm doesn't require a separate indexing step)

Align*Pass* #Scripts to align RNAseq reads to genome using STAR in 2-pass mode. Pass 1 - align. Pass 2 - align with splice junction databases.

		Align*Pass* #Scripts to align RNAseq reads to genome using STAR in 2-pass mode

		NGMLoopFelixPopHiC11-1-18 #Script to align DNAseq reads to genome using NextGenMapping


	#### Picard_process_alignments

		Readgroups_duplicates_sort_\* #PicardTools scripts to post-process alignments

	#### SNP_calls
		
		MPileup\* and Call\* #Scripts to generate genotype calls from alignments using bcftools 1.9-67-g626e46b

	#### Filter\
	
		Filter\* #Scripts to filter genotype calls by quality and depth
		
	#### Convert to 012
		
		Sample scripts to convert data to 012 format using VCFtools

#### Datasets 

ComparativeTranscriptomeLinkageMappingTXF2s - TX F2 mapping population, RNAseq
ComparativeTranscriptomeLinkageMappingNCF2s - NC F2 mapping population, RNAseq
FelixPop - population data from Beaudry FEG, Barrett SCH, Wright SI. 2019. Ancestral and neo‐sex chromosomes contribute to population divergence in a dioecious plant. Evolution. 54:180, RNASeq
JoshPop - population data from Hough J, Hollister JD, Wang W, Barrett SCH, Wright SI. 2014. Genetic degeneration of old and young Y chromosomes in the flowering plant Rumex hastatulus. Proc Natl Acad Sci. 111(21):7713–7718, RNASeq
JoshCross - pedigree cross data from Hough J, Hollister JD, Wang W, Barrett SCH, Wright SI. 2014. Genetic degeneration of old and young Y chromosomes in the flowering plant Rumex hastatulus. Proc Natl Acad Sci. 111(21):7713–7718, RNASeq
FelixRADSEQPop, popradseq - population data from Beaudry FEG, Barrett SCH, Wright SI. 2017. Genomic loss and silencing on the Y chromosomes of Rumex. Genome Biol Evol. 9(12):3345–3355m DNAseq
