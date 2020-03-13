
print("hello")

import re


agp_file_path = input("Specify path to agp file: ")
marker_file_path = input("Specify path to marker file: ")
outfile_path = input("Specify path to output file: ")
#marker_file = ("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/NA_RM_mapped_and_colocated_markers_TX.csv")

##New map
#agp_file_path = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp
#marker file = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/Final_markers_for_conversion.csv
#output file E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/orientation_aware_AGP_converted.csv

#agp_file_path = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp
#marker file = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/Final_markers_for_conversion.csv
#output file D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/newAGP_converted.csv

##old map
#agp_file_path = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/oldCHRR_genome.agp
#marker file = D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/NA_RM_mapped_and_colocated_markers_TX.csv
#output file D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/AGP_converted_old_map.csv


agp={}

#print(marker_file)
linecounter = 0
with open(agp_file_path, "r") as f:
    for line in f:
        if line[0] != "#":
            print(line)
            line = line.strip()
            values = re.split('\t', line)
            #print(values)
            #print(values[5])
            #print(values[5][0])
            if values[5][0] == "S":
                chrrLG = values[0]
                print(chrrLG)
                LG_pos = (int(values[1]), int(values[2]))
                print(LG_pos)
                scaff_pos = (int(values[6]), int(values[7]))
                print(scaff_pos)
                orientation=values[8]
                print(orientation)
                scaffold=int(re.split('\_', values[5])[1])
                print(scaffold)
                #agp[scaffold].append(scaff_pos, LG_pos, orientation)
                #Stuck on making list of tuples that is not class tuple
                if scaffold not in agp.keys():
                    #agp[scaffold]={}
                    #agp[scaffold]=list(scaff_pos, LG_pos, orientation)
                    agp[scaffold]=[]
                    print("key")
                    print(agp)
                    print(agp[scaffold])
                    agp[scaffold].append((scaff_pos, LG_pos, orientation, chrrLG))
                    #agp[scaffold][scaff_pos]=(LG_pos, orientation)
                else:
                    agp[scaffold].append((scaff_pos, LG_pos, orientation, chrrLG))
        linecounter += 1
        print("line count ", linecounter)

print(agp)
values=[]

#marker_file = open(agp_file_path, "r")
#marker_file = ("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/short_NA_RM_mapped_and_colocated_markers_TX.csv")
#marker_file = ("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/NA_RM_mapped_and_colocated_markers_TX.csv")
#agp={}
#print(marker_file)
linecounter = 0
outfile=[]
with open(marker_file_path, "r") as f:
    for line in f:
        if linecounter>0:
        #if line[0] != "#":
            #print(line)
            line = line.strip()
            values = re.split(',', line)
            #print(values)
            scaffold=int(values[1])
            #print(scaffold)
            #print(type(scaffold))
            position = int(values[2])
            ASMAP_LG = int(values[3])
            CM = float(values[4])
            #print(CM)
            print(linecounter)
            #print(values)
            if scaffold in agp.keys():
                print("FOUND")
                for i in agp[scaffold]:
                    print(i)

                    #print(i[1][0])
                    #print(i[0][0])
                    #print(i[0][1])
                    #print(type(i[0][1]))
                    #print(position)
                    #print(type(position))
                  #  print(type(position))
                   # print((i[0][0] < position ))#< (i[1][1])))
                   # print((position < i[0][1]))
                   # print(i[0][0] < position < i[0][1])
               #print("RANGED")
                   # new_position="NA"
                    if (i[0][0] < position < i[0][1])==True:
                        print("YIPPEE")
                        print(position)
                        print(i[2])
                        LG=i[3]

                        #new_position = (((i[0][1]) - position) + (i[1][0]))
                       # print((i[0][1]), (i[1][0]))
                       # print("new position ", new_position)
                       # agp[scaffold].append((scaff_pos, LG_pos, orientation, chrrLG))

                        if (i[2]=='+'):
 #                           #  = position on scaffold fragment - scaffold fragment start + LG section start
                            print("Forward")
                            #                        new_position=(position-i[0][0]+i[1][0])
                            new_position=(((position-i[0][0]))+(i[1][0]))
 #                           print((i[0][1]),(i[1][0]))
 #                           print("new position ", new_position)

                        elif (i[2] == '-'):
 #                           # = scaffold fragment end - position on scaffold + LG section start
 #                           print("Reverse")
                            new_position=((i[0][1]-position)+i[1][0])
 #                           print((i[0][1]),(i[1][0]))
 #                           print("new position ", new_position)
  #                      else:
  #                          new_position="NA"
            values.append(new_position)
            values.append(LG)
            outfile.append(values)
        linecounter += 1
#print(outfile)
 
linecounter=0

#f = open("AGP_converted.csv", "w+")
f = open(outfile_path, "w+")

f.write("Marker_number, Scaffold, Position_on_scaffold, LG_in_ASMAP, CM_position, Position_on_LG, LG,\n")
#6830,17266,33348559,5,1.818466791,5650015,L.1,
print(outfile[0])
for i in range(0,len(outfile)):
#    print((outfile[i]))
#    print(str((outfile[i])))
    for j in (outfile[i]):
        f.write(str(j))
        f.write(",")
        print(str(j))
    f.write("\n")
f.close()
