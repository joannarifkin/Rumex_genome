library(dplyr)
library(tidyverse)
scipen=999
setwd("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses")

#Import windowed files, output from Windowed_analysis_no_conversion_ohta_version

#windowed_table_file<-"TX_full_map_5000000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"TX_full_map_1000000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"TX_full_map_500000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"TX_full_map_100000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"TX_full_map_50000_bp_recombination_rates_2020-05-27.csv"


#windowed_table_file<-"NC_full_map_5000000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"NC_full_map_1000000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"NC_full_map_500000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"NC_full_map_100000_bp_recombination_rates_2020-05-27.csv"
#windowed_table_file<-"NC_full_map_50000_bp_recombination_rates_2020-05-27.csv"

#From here to end, run for each block of windowed import files (NC and TX)
windowed_table<-read.csv(windowed_table_file, stringsAsFactors = F)

header<-(strsplit(windowed_table_file, "_"))
window_size<-header[[1]][4]
cytotype<-header[[1]][1]

#25 MB windows
MB25_summary_windowed_table<-windowed_table %>% group_by(LG) %>% mutate(LG_block=Position_on_LG %/%25000000)%>%group_by(LG,LG_block)%>%select_if(., is.numeric) %>% summarize_all(funs(mean(., na.rm=T)))%>%mutate(block_size="25MB")



#50 MB windows
MB50_summary_windowed_table<-windowed_table %>% group_by(LG) %>% mutate(LG_block=Position_on_LG %/%50000000)%>%group_by(LG,LG_block)%>%select_if(., is.numeric) %>% summarize_all(funs(mean(., na.rm=T)))%>%mutate(block_size="50MB")
View(MB50_summary_windowed_table)

#25 CM windows
CM25_summary_windowed_table<-windowed_table %>% group_by(LG) %>% mutate(LG_block=CM_position %/%25)%>%group_by(LG,LG_block)%>%select_if(., is.numeric) %>% summarize_all(funs(mean(., na.rm=T)))%>%mutate(block_size="25CM")


#50 CM windows
CM50_summary_windowed_table<-windowed_table %>% group_by(LG) %>% mutate(LG_block=CM_position %/%50)%>%group_by(LG,LG_block)%>%select_if(., is.numeric) %>% summarize_all(funs(mean(., na.rm=T)))%>%mutate(block_size="50CM")

Windowed_tables<-rbind(MB25_summary_windowed_table, MB50_summary_windowed_table, CM25_summary_windowed_table, CM50_summary_windowed_table)

Windowed_tables<-Windowed_tables %>% mutate(cytotype=cytotype) %>% mutate(window_size=window_size)
#View(All_windowed_tables)

#All_windowed_tables<-Windowed_tables ##Use this for the first file
All_windowed_tables<-rbind(All_windowed_tables, Windowed_tables) #Use this for subsequent files

#End of block to repeat for all five NC or TX files

All_windowed_tables<-All_windowed_tables%>%mutate(coefficient_MB=coefficient*1000000)
All_windowed_tables<-All_windowed_tables%>%mutate(LG_final=str_replace_all(LG, c("L.10"="Sex",
"L.7"="A1","L.8"="A2","L.3"="A4","L.5"="A3")))
  

write.csv(All_windowed_tables, paste(as.character(cytotype), "_window_blocks_summary.csv", sep=""), quote=F)



