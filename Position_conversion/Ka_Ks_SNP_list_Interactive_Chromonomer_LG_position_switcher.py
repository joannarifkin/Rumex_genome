

print("hello")


import re
from operator import itemgetter



agp_file_path = input("Specify path to agp file: ")
snp_file_path = input("Specify path to SNP file: ")
outfile_path = input("Specify path to output file: ")


#agp_file_path = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/FinalCHRR_genome.agp
#snp_file_path = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/Divergence_data/Ka_Ks_positions_only.csv
#out_file_path = E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/Divergence_data/Ka_Ks_positions_only_converted.csv


agp={}
#print(marker_file)
linecounter = 0
with open(agp_file_path, "r") as f:
    for line in f:
        if line[0] != "#":
            #print(line)
            line = line.strip()
            values = re.split('\t', line)
            print(values)
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
               # newscaff = (scaffparts[0] + "_" + scaffparts[1] + ";" + scaffparts[2] + "=" + scaffparts[3])
                newscaff = int(scaffparts[1])
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
print(header)
with open(snp_file_path, "r") as f:
    for line in f:
        newposition=0
        print(line)
        #if line == 0:
        #    header=line
        if linecounter > 0:
            line = line.strip()
            values = re.split(',', line)
            print(values)
            outline = []
            snp_scaffold=int(values[0])
            position_1 = int(values[1])
            position_2 = int(values[2])
            print("scaffold", snp_scaffold)
            print("start", position_1)
            print("end", position_2)
            #print(content)
            if snp_scaffold in agp.keys():
                print(snp_scaffold)
                print("found")
                for i in agp[snp_scaffold]:
                    #print(i)
                    if (i[0][0] < position_1 < i[0][1])==True:
                        LG = i[3]
                        print(LG)
                        if (i[2]=='+'):
                            new_position_1=(((position_1-i[0][0]))+(i[1][0]))
                        elif (i[2] == '-'):
                            new_position_1=((i[0][1]-position_1)+i[1][0])
                        print(new_position_1)
                    if (i[0][0] < position_2 < i[0][1]) == True:
                        LG = i[3]
                        print(LG)
                        if (i[2] == '+'):
                            new_position_2 = (((position_2 - i[0][0])) + (i[1][0]))
                        elif (i[2] == '-'):
                            new_position_2 = ((i[0][1] - position_2) + i[1][0])
                        print(new_position_2)
            else:
                LG = "NA" #Because NA can't be sorted
                new_position=0 #Because NA can't be sorted
         #print(linecounter)
        #linecounter += 1
            print("in scaffold", snp_scaffold)
            outline.append(snp_scaffold)
            print("in start", position_1)
            outline.append(position_1)
            print("in end", position_2)
            outline.append(position_2)
            print("out", LG)
            outline.append(LG)
            print("out start", new_position_1)
            outline.append(new_position_1)
            print("out end", new_position_2)
            outline.append(new_position_2)
            outfile.append(outline)
        else:
            line = line.strip()
            values = re.split(',', line)
            print(values)
            #line = line.strip()
            header.append(values[0])
            header.append(values[1])
            header.append(values[2])
            header.append("LG")
            header.append("LG_position_start")
            header.append("LG_position_end")
        linecounter += 1
# print(header)
print(len(outfile))
print(outfile)
print(type(outfile[1][3]))

outfile_sorted=sorted(outfile, key=itemgetter(3,4))
#sorted(L, key=itemgetter(2))
print(outfile_sorted)
linecounter=0
#f = open("AGP_converted.csv", "w+")
f = open(outfile_path, "w+")
#f.write("Marker_number, Scaffold, Position_on_scaffold, LG_in_ASMAP, CM_position, Position_on_LG, LG,\n")
#6830,17266,33348559,5,1.818466791,5650015,L.1,
#print(outfile[0])

for i in range (0, len(header)):
    f.write(header[i])
   # print(header[i])
    f.write(",")
f.write("\n")
#f.write(header)
#print(header)
for i in range(0,len(outfile_sorted)):
#    print((outfile[i]))
#    print(str((outfile[i])))
    for j in range(0,(len(outfile_sorted[i])-1)):
        f.write(str(outfile_sorted[i][j]))
        f.write(",")
      #  print(str(j))
    f.write(str(outfile_sorted[i][-1]))
    f.write("\n")
#print(str(outfile[-1]))
f.close()
