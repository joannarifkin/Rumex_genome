/ohta/joanna.rifkin/bcftools/bcftools mpileup -d 3000 -Ou -f /ohta/joanna.rifkin/Genomes/Hi-C/hastate_28Sep2018_nbKXS.fasta /ohta/joanna.rifkin/Alignments/STAR/JoshHastatulusCross/Aligned_To_HiC_Assembly/Picard/AnalysisReady_*.bam | /ohta/joanna.rifkin/bcftools/bcftools call -vmO z -f gq -o /ohta/joanna.rifkin/HiCSNPData/JoshCross/raw_snps_11-9-18.vcf.gz


 -a, --annotate LIST

    Comma-separated list of FORMAT and INFO tags to output. (case-insensitive, the "FORMAT/" prefix is optional, and use "?" to list available annotations on the command line) [null]:

    *FORMAT/AD* .. Allelic depth (Number=R,Type=Integer)
    *FORMAT/ADF* .. Allelic depths on the forward strand (Number=R,Type=Integer)
    *FORMAT/ADR* .. Allelic depths on the reverse strand (Number=R,Type=Integer)
    *FORMAT/DP* .. Number of high-quality bases (Number=1,Type=Integer)
    *FORMAT/SP* .. Phred-scaled strand bias P-value (Number=1,Type=Integer)

    *INFO/AD* .. Total allelic depth (Number=R,Type=Integer)
    *INFO/ADF* .. Total allelic depths on the forward strand (Number=R,Type=Integer)
    *INFO/ADR* .. Total allelic depths on the reverse strand (Number=R,Type=Integer)

    *FORMAT/DV* .. Deprecated in favor of FORMAT/AD;
            Number of high-quality non-reference bases, (Number=1,Type=Integer)
    *FORMAT/DP4* .. Deprecated in favor of FORMAT/ADF and FORMAT/ADR;
            Number of high-quality ref-forward, ref-reverse,
            alt-forward and alt-reverse bases (Number=4,Type=Integer)
    *FORMAT/DPR* .. Deprecated in favor of FORMAT/AD;
            Number of high-quality bases for each observed allele (Number=R,Type=Integer)
    *INFO/DPR* .. Deprecated in favor of INFO/AD;
            Number of high-quality bases for each observed allele (Number=R,Type=Integer)