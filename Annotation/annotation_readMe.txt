This annotation is based on BRAKER and gene prediction from Arabidopsis. Prediction and filtering is informed with diverse Rumex hastatulus resources.

Steps
1. A repeat library is made and the assembly is masked for repeats.
2. RNAseq is aligned to the assembly using Star
3. BRAKER uses the masked assembly, the gene library from Arabidopsis and RNAseq to predict genes.
4. BRAKER output had a very high count of predicted genes, likely due to high TE density. We therefore filter (with R script) the predicted genes by:
	4a. we mapped the RNAseq data back to the predicted transcriptome to only include genes that are expressed. 
	4b. we mapped whole genome sequences back to the predicted transcriptome to exclude genes with more coverage than would be predicted for single copy genes.
5. The .gff was than filetered for genes that passed the above (4) threshold