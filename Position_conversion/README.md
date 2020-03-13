# Position_conversion
Fixed_Marker_fasta_generator.py #Python script to extract SNP positions from the genome as 100-bp fragments for alignment to produce Chromonomer input 
Terminal_input.txt #Inputs for running Chromonomer on the command line 
"LGPositionSwitcher" scripts #All these scripts follow the same basic format: they read a .agp file from Chromonomer into a dictionary with scaffolds as keys and linkage map positions as entries. They then read in an input file (.gff, .vcf, .txt, .csv, etc.) with positions based on the scaffolds and output new positions in the linkage groups. 