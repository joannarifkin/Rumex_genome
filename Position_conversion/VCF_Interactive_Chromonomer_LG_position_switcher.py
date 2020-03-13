

print("hello")


import re


agp_file_path = input("Specify path to agp file: ")
vcf_file_path = input("Specify path to vcf file: ")
outfile_path = input("Specify path to output file: ")


#agp_file_path = "D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp"
#vcf_file_path = "D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/miniVCF.vcf"
#outfile_path = "testout.vcf"


#marker_file = ("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/NA_RM_mapped_and_colocated_markers_TX.csv")
#VCF "D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/miniVCF.vcf"
#VCF "E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/test.vcf"
##New map
#agp_file_path = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp
#marker file = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/Final_markers_for_conversion.csv
#output file E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/VCF_AGP_converted.csv
#output file D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/VCF_AGP_converted.csv

#agp_file_path = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp
#marker file = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/Final_markers_for_conversion.csv
#output file D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/newAGP_converted.csv

agp={}
#print(marker_file)
linecounter = 0
with open(agp_file_path, "r") as f:
    for line in f:
        if line[0] != "#":
            #print(line)
            line = line.strip()
            values = re.split('\t', line)
            #print(values)
            #print(values[5])
            #print(values[5][0])
            if values[5][0] == "S":
                chrrLG = values[0]
                #print(chrrLG)
                LG_pos = (int(values[1]), int(values[2]))
                #print(LG_pos)
                scaff_pos = (int(values[6]), int(values[7]))
                #print(scaff_pos)
                orientation=values[8]
                #print(orientation)
                #Need to fix scaffold so there is a semicolon rather than an underscore
                scaffold=values[5]
                scaffparts = re.split('\_', scaffold)
                #print(scaffold)
                #print(scaffparts)
                newscaff = (scaffparts[0] + "_" + scaffparts[1] + ";" + scaffparts[2] + "=" + scaffparts[3])
                #print(newscaff)
                #agp[scaffold].append(scaff_pos, LG_pos, orientation)
                #Stuck on making list of tuples that is not class tuple
                if newscaff not in agp.keys():
                    #agp[scaffold]={}
                    #agp[scaffold]=list(scaff_pos, LG_pos, orientation)
                    agp[newscaff]=[]
                    #print("key")
                    #print(agp)
                    #print(agp[newscaff])
                    agp[newscaff].append((scaff_pos, LG_pos, orientation, chrrLG))
                    #agp[scaffold][scaff_pos]=(LG_pos, orientation)
                else:
                    agp[newscaff].append((scaff_pos, LG_pos, orientation, chrrLG))

        linecounter += 1
        print("line count ", linecounter)
#print(header)
#
#print(agp)

values=[]
linecounter = 0
outfile=[]
header=[]
with open(vcf_file_path, "r") as f:
    for line in f:
        newposition=0
        #print(line)
        if line[0] != "#":
            outline = []
            line = line.strip()
            values = re.split('\t', line)
            scaffold_vcf=values[0]
            position = int(values[1])
            content=values[2:]
            print(scaffold_vcf)
            print(position)
            #print(content)
            if scaffold_vcf in agp.keys():
                print(scaffold_vcf)
                print("found")
                for i in agp[scaffold_vcf]:
                    #print(i)
                    if (i[0][0] < position < i[0][1])==True:
                        LG = i[3]
                        print(LG)
                        if (i[2]=='+'):
                            new_position=(((position-i[0][0]))+(i[1][0]))
                        elif (i[2] == '-'):
                            new_position=((i[0][1]-position)+i[1][0])
                        print(new_position)
                print(scaffold_vcf, "out")
                outline.append(scaffold_vcf)
                print(position, "out")
                outline.append(position)
                print(LG, "out")
                outline.append(LG)
                print(new_position, "out")
                outline.append(new_position)
                for i in range(0,len(content)):
                    #print(content[i])
                    outline.append(content[i])
                outfile.append(outline)
        else:
            line = line.strip()
            header.append(line)
        linecounter += 1
#print(outfile)
 #   print(content)
  #  print(outfile)
# print(header)
print(len(outfile))

linecounter=0
#f = open("AGP_converted.csv", "w+")
f = open(outfile_path, "w+")
#f.write("Marker_number, Scaffold, Position_on_scaffold, LG_in_ASMAP, CM_position, Position_on_LG, LG,\n")
#6830,17266,33348559,5,1.818466791,5650015,L.1,
#print(outfile[0])
for i in range (0, len(header)):
    f.write(header[i])
   # print(header[i])
    f.write("\n")
#f.write(header)
#print(header)
for i in range(0,len(outfile)):
#    print((outfile[i]))
#    print(str((outfile[i])))
    for j in range(0,(len(outfile[i])-1)):
        f.write(str(outfile[i][j]))
        f.write("\t")
      #  print(str(j))
    f.write(str(outfile[i][-1]))
    f.write("\n")
#print(str(outfile[-1]))
f.close()
