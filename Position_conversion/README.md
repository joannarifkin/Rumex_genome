# Position_conversion
Fixed_Marker_fasta_generator.py #Python script to extract SNP positions from the genome as 100-bp fragments for alignment to produce Chromonomer input 
Terminal_input.txt #Inputs for running Chromonomer on the command line 

# Position conversion scripts

#All these scripts follow the same basic format: they read a .agp file from Chromonomer into a dictionary with scaffolds as keys and linkage map positions as entries. They then read in an input file (.gff, .vcf, .txt, .csv, etc.) with positions based on the scaffolds and output new positions in the linkage groups. 

VCF_Interactive_Chromonomer_LG_position_switcher.py #Converts a VCF (used to convert VCFs for LDHat and GWAS analyses)
gff_Interactive_Chromonomer_LG_position_switcher.py #Converts a GFF (used for windowing annotation data)
Fixed_Marker_fasta_generator.py #Generates 100-bp markers as input for Chromonomer
Ka_Ks_SNP_list_Interactive_Chromonomer_LG_position_switcher.py #Converts list of positions with Ka/Ks values from Hough J, Hollister JD, Wang W, Barrett SCH, Wright SI. 2014. Genetic degeneration of old and young Y chromosomes in the flowering plant Rumex hastatulus. Proc Natl Acad Sci. 111(21):7713–7718. 
Map_file_SNP_list_Interactive_Chromonomer_LG_position_switcher.py #Converts linkage map markers in a .csv file (used for several analyses)