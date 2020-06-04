library(dplyr)
library(stringr)
library(tidyr)
setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits")
setwd("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits")

###### Add colocated markers to map ---------------------

library(dplyr)
library(tidyr)
full_map<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/Final_F2_style_TX_map_biggest_100_scaffolds_plus_sex_linked_10-21.csv", stringsAsFactors = F) #Import map exported from ASMAP
full_map<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/Final_F2_style_TX_map_biggest_100_scaffolds_plus_sex_linked_10-21.csv", stringsAsFactors = F) #Import map exported from ASMAP
full_map<-full_map[,1:3] #Remove individual genotypes
full_colocated<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/Colocated_markers_biggest_100_plus_sex_linked.csv", stringsAsFactors = F) #Import colocated marker bins exported from ASMAP
full_colocated<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/Colocated_markers_biggest_100_plus_sex_linked.csv", stringsAsFactors = F) #Import colocated marker bins exported from ASMAP



full_map<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/NC_transcriptome/F2_style_NC_map_clean_biggest100_5-28.csv", stringsAsFactors = F) #Import map exported from ASMAP
#full_map<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/PleaseLetThisBeItNC_map_biggest_100_scaffolds_sex_linked_split_11-8.csv", stringsAsFactors = F) #Import map exported from ASMAP
full_map<-full_map[,1:3] #Remove individual genotypes
full_colocated<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/NC_transcriptome/Colocated_markers_biggest_100.csv", stringsAsFactors = F) #Import colocated marker bins exported from ASMAP
#full_colocated<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_construction/NCColocated_markers_biggest_100_plus_sex_linked.csv", stringsAsFactors = F) #Import colocated marker bins exported from ASMAP





#head(full_colocated)
#head(full_map)
#length(intersect(full_map[,1], full_colocated$mark)) #Identify number of markers in map with associated bins of colocated markers
#View(intersect(full_map[,1], full_colocated$mark))

#head(full_map)
positions_colocated<-full_join(full_colocated, full_map,  by=c("mark"="Genotype"), fill=T) #Join list of colocated markers to map
View(subset(positions_colocated, positions_colocated$bins=="239"))

positions_colocated_filled<-positions_colocated%>%group_by(bins)%>%fill(X.y,X.1) #Fill centimorgan and linkage group info for each bin of colocated markers

positions_colocated_filled<-subset(positions_colocated_filled,is.na(positions_colocated_filled$X.1)==F) #Discard markers that were not placed in the map
head(positions_colocated_filled)
positions_colocated_filled<-positions_colocated_filled[,-3]
colnames(positions_colocated_filled)<-c("index","bin","marker","LG","CM")
write.csv(positions_colocated_filled,"Old_Map_NC_Final_markers_positions_colocated_filled.csv", row.names=T, quote = F)

#################### Move small scaffolds from small linkage groups to parallel full linkage groups TX #######

positions_colocated_filled<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Final_markers_positions_colocated_filled.csv", stringsAsFactors = F, row.names=1)
positions_colocated_filled<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Final_markers_positions_colocated_filled.csv", stringsAsFactors = F, row.names=1)
#Move markers from unique scaffolds on LG4 to LG10

positions_colocated_filled<-data.frame(positions_colocated_filled)
colnames(positions_colocated_filled)
str(positions_colocated_filled)
positions_colocated_filled<-separate(positions_colocated_filled, marker, into=c(NA, "scaffold", "position"), se=("\\_"), remove = F) %>% separate(., scaffold, into=c("scaffold", NA, NA), se=c("\\-"), remove = F)


CM_4_10_difference<-(subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_17040-HRSCAF-21160_3444680")$CM)-
  (subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_17040-HRSCAF-21160_4391861")$CM)

positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      LG == "L.4" ~ (CM+CM_4_10_difference ),
      TRUE ~ CM)) %>%
  mutate (
    LG = case_when (
      LG == "L.4" ~ "L.10",
      TRUE ~ LG))
 
CM_15769<- (subset (positions_colocated_filled, positions_colocated_filled$marker=="ScnbKXS_17265-HRSCAF-21557_12838208"))$CM
positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      marker == "ScnbKXS_15769-HRSCAF-18489_111981" ~ (CM_15769),
      TRUE ~ CM)) %>%
  mutate (
    LG = case_when (
      marker == "ScnbKXS_15769-HRSCAF-18489_111981" ~ "L.5",
      TRUE ~ LG))

#Fix LG 11 unique scaffold positions
CM_292<-(subset (positions_colocated_filled, positions_colocated_filled$marker=="ScnbKXS_15927-HRSCAF-18817_297474"))$CM
#CM_3947<-(subset (positions_colocated_filled, positions_colocated_filled$marker=="ScnbKXS_15927-HRSCAF-18817_297474"))$CM #These two scaffolds are at the same genetic position
CM_14219<-(subset (positions_colocated_filled, positions_colocated_filled$marker=="ScnbKXS_17162-HRSCAF-21394_26890388"))$CM

positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      (scaffold == 292 | scaffold == 3947) ~ (CM_292),
      TRUE ~ CM)) %>%
  mutate (
    LG = case_when (
      (scaffold == 292 | scaffold == 3947) ~ "L.7",
      TRUE ~ LG))


positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      (scaffold == 14219) ~ (CM_14219),
      TRUE ~ CM)) %>%
  mutate (
    LG = case_when (
      (scaffold == 14219) ~ "L.7",
      TRUE ~ LG))

#View(test)

positions_colocated_filled<-filter(positions_colocated_filled, !(LG %in% c("L.1", "L.2","L.4","L.6","L.9","L.11", "L.12")))

write.csv(positions_colocated_filled, "Final_edited_markers_positions_colocated_filled.csv", quote = F, row.names = F)
  


#################### Move small scaffolds from small linkage groups to parallel full linkage groups NC #######

positions_colocated_filled<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Old_Map_NC_Final_markers_positions_colocated_filled.csv", stringsAsFactors = F, row.names=1)
#positions_colocated_filled<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/NC_Final_markers_positions_colocated_filled.csv", stringsAsFactors = F, row.names=1)



positions_colocated_filled<-data.frame(positions_colocated_filled)
colnames(positions_colocated_filled)
str(positions_colocated_filled)
positions_colocated_filled<-separate(positions_colocated_filled, marker, into=c(NA, "scaffold", "position"), se=("\\_"), remove = F) %>% separate(., scaffold, into=c("scaffold", NA, NA), se=c("\\-"), remove = F)




#17271, 17272, 4250, 6377 - 17040	4008573	L.10

positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      LG == "L.11" & (scaffold == 4250 | scaffold == 6377 | scaffold == 17271 | scaffold == 17272)      ~ (subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_17040-HRSCAF-21160_4008573")$CM),
      TRUE ~ CM))  %>%
  mutate (
    LG = case_when (
      LG == "L.11" & (scaffold == 4250 | scaffold == 6377 | scaffold == 17271 | scaffold == 17272)   ~ "L.10",
      TRUE ~ LG))

#7031 - 17265	18616492 L.10


#8361 - 11619	3218474	L.10 - ((11619	3213180	L.5)-(8361	116661	L.5))


positions_colocated_filled<-positions_colocated_filled %>%
  mutate (
    CM = case_when (
      LG == "L.5" & (scaffold == 8361) ~ (subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_11619-HRSCAF-13537_3218474")$CM - (((subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_11619-HRSCAF-13537_3213180")$CM) - ((subset(positions_colocated_filled, positions_colocated_filled$marker == "ScnbKXS_8361-HRSCAF-9758_116661")$CM))))),
      TRUE ~ CM))  %>%
  mutate (
    LG = case_when (
      LG == "L.5" & (scaffold == 8361) ~ "L.10",
      TRUE ~ LG))


View(positions_colocated_filled)
positions_colocated_filled<-filter(positions_colocated_filled, (LG %in% c("L.1", "L.7","L.9","L.10")))

write.csv(positions_colocated_filled, "Old_Map_NC_Final_edited_markers_positions_colocated_filled.csv", quote = F, row.names = F)




################ Reformat map and marker files for chromonomer ####################


setwd("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Chromonomer_final")
head(positions_colocated_filled)

Chromonomer_input_map<-positions_colocated_filled[,c(6,3,7)] #LG, marker, position
write.table(Chromonomer_input_map, "TX_chromonomer_input_map_final.txt", quote=F, row.names = F,col.names = F, sep = "/t", eol = "\n")

Chromonomer_marker_names<-positions_colocated_filled[,c(3)]
write.table(Chromonomer_marker_names, "TX_chromonomer_markers_final.txt", quote=F, row.names = F,col.names = F, eol = "\n")


Chromonomer_input_map<-positions_colocated_filled[,c(6,3,7)] #LG, marker, position
write.table(Chromonomer_input_map, "Old_Map_NC_chromonomer_input_map_final.txt", quote=F, row.names = F,col.names = F, sep = "\t", eol = "\n")

Chromonomer_marker_names<-positions_colocated_filled[,c(3)]
write.table(Chromonomer_marker_names, "Old_Map_NC_chromonomer_markers_final.txt", quote=F, row.names = F,col.names = F, eol = "\n")


############## Format markers for chromonomer AGP conversion script ###################

setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Chromonomer_final")
library(tidyr);library(dplyr)

positions_colocated_filled<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Final_edited_markers_positions_colocated_filled.csv", header=T, stringsAsFactors = F, row.names = F)
positions_colocated_filled<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Final_edited_markers_positions_colocated_filled.csv", header=T, stringsAsFactors = F)
colnames(positions_colocated_filled)
#conversion_markers<-positions_colocated_filled %>% separate(marker, int0)


#positions_colocated_filled<-separate(positions_colocated_filled, marker, into=c(NA, "scaffold", "position"), se=("//_"), remove = F) %>% separate(., scaffold, into=c("scaffold", NA, NA), se=c("//-"), remove = F)


markers_for_conversion<-positions_colocated_filled[,c(5,4,6, 7)]
colnames(markers_for_conversion)
markers_for_conversion<-separate(markers_for_conversion, LG, into=c(NA, "LG"), sep="\\.")

#write.csv(markers_for_conversion, "TXFinal_markers_for_conversion.csv", quote = F)
#write.csv(markers_for_conversion, "Old_Map_NCFinal_markers_for_conversion.csv", quote = F)


colnames(full_map)
lengths<-full_map %>% group_by(X) %>% summarize(lglength=max(X.1))

View(lengths)
write.table(lengths, "lengths")

