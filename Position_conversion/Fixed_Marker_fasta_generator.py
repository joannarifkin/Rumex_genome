

print("hello")
import re



marker_file_path = input("Specify path to marker file: ")
genome_file_path = input("Specify path to genome file: ")
outfile_path = input("Specify path to output file: ")
#("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/LepMap/TX_transcriptome/practice_post_file_modified", "w")




genome={}
genome_file = open(genome_file_path, "r")
#print(marker_file)
with open(genome_file_path, "r") as f:
    linecounter=0
    for line in f:
        print(line[0])
        if line[0]==">":
            line=line.strip()
            values = re.split('_|=|;', line)
            #print(values)
            key=values[1]
            genome[key] = (next(f).strip()).upper()
            #print(key)
            #print(genome[key])
print(len(genome.keys()))


#print(genome[3])
print(marker_file_path)
marker_file = open(marker_file_path, "r")
#print(marker_file)

markers_out={}
for line in marker_file:
    print(line)
    line = line.strip()
    print(line)
    values=[]
    values = re.split('_|-', line)
    #    values = line.split("\_")
    print(values)
    scaffold=values[1]
    print(scaffold)
    position=values[4]
    print(position)
    print(type(position))
    #print(genome[scaffold])
    print(genome[scaffold][int(position)-1])
    print(genome[scaffold][int(position):int(position) + 100])
    key = (">" + values[0] + "_" + values[1] + "-" + values[2] + "-" + values[3] + "_" + values[4])
    print(key)
    markers_out[key]=genome[scaffold][int(position):int(position)+100]
    print(markers_out[key])

print(markers_out)


with open(outfile_path, 'w') as f:
    for key in markers_out:
        print(key, file=f)
        print(markers_out[key], file=f)
