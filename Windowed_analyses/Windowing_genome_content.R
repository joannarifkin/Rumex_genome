library(dplyr)
library(tidyr)
#install.packages("ff")
library(ff)
options(scipen=999)

#setwd("/ohta/joanna.rifkin/HiCSNPData/Windowed_analyses/Genome_content/")
getwd()



################# Import data ###########

#All_genes<-read.table("/ohta/felix.beaudry/beyond/ann2/goodGenes.chrom.gff", stringsAsFactors = F)



#/ohta/felix.beaudry/beyond/ann2/goodGenes.chrom.gff
head(All_genes)
colnames(All_genes)[1:5]<-c("LG","Source", "Class", "LG_start","LG_end")
All_genes<-mutate(All_genes, LG_lower_start=pmin(LG_start,LG_end))
All_genes<-mutate(All_genes, featuresize=abs(LG_start-LG_end))
All_genes<-subset(All_genes, All_genes$Class == "gene")
head(All_genes)

All_repeats<-read.csv("/ohta/joanna.rifkin/HiCSNPData/Windowed_analyses/Genome_content/JR_LG_repeats.gff", stringsAsFactors = F)
All_repeats<-mutate(All_repeats, LG_lower_start=pmin(LG_start,LG_end))
All_repeats<-mutate(All_repeats, featuresize=abs(LG_start-LG_end))

head(All_repeats)
unique(All_repeats$transcript_status )
repeat_distribution<-All_repeats%>%group_by(LG,transcript_status)%>% tally
write.csv(repeat_distribution, "Repeats_by_class.csv")
#View(All_genes$start)
#hist(All_genes_BUSCO$start)
#hist(All_repeats$start)

######### Join data #####

All_features<-rbind(All_repeats, All_genes)
head(All_features)
All_features<-arrange(All_features, LG, LG_start)
#View(All_features)

All_features<-distinct(All_features)

write.csv(All_features, "All_features.csv")

unique(All_features$LG)
All_features_L10<-subset(All_features, All_features$LG=="L.10")
c = outer(All_features$LG_start, All_features$LG_end, ">")
d = outer(All_features$LG_end, All_features$LG_start, "<")



for (i in (unique(All_features$LG)))
{
  print(i)
  All_features_sub
  
  
}

gc()
c = ff(outer(All_features$LG_start, All_features$LG_end, ">"))
d = outer(All_features$LG_end, All_features$LG_start, "<")

All_features <- All_features %>% 
  mutate(Overlap = apply(c & d, 1, sum) > 1 
  )

######### Window data #####
### 10000
All_genes_windowed_10K <- All_genes %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%10000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(All_genes_windowed_10K)
All_genes_windowed_10K<-All_genes_windowed_10K[,c(1:2,6:8,11)]
All_genes_windowed_10K<-All_genes_windowed_10K %>% mutate(proportion_of_window=featuresize_sum/10000)%>%mutate(window_start=(position_window*10000+1))
write.csv(All_genes_windowed_10K, "6_2_All_genes_windowed_10K_with_feature_size.csv", row.names=F, quote=F)

### 50000
All_genes_windowed_50K <- All_genes %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%50000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(All_genes_windowed_50K)
All_genes_windowed_50K<-All_genes_windowed_50K[,c(1:2,6:8,11)]
All_genes_windowed_50K<-All_genes_windowed_50K %>% mutate(proportion_of_window=featuresize_sum/50000)%>%mutate(window_start=(position_window*50000+1))
write.csv(All_genes_windowed_50K, "6_2_All_genes_windowed_50K_with_feature_size.csv", row.names=F, quote=F)

### 100000
All_genes_windowed_100K <- All_genes %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%100000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(All_genes_windowed_100K)
All_genes_windowed_100K<-All_genes_windowed_100K[,c(1:2,6:8,11)]
All_genes_windowed_100K<-All_genes_windowed_100K %>% mutate(proportion_of_window=featuresize_sum/100000)%>%mutate(window_start=(position_window*100000+1))
write.csv(All_genes_windowed_100K, "6_2_All_genes_windowed_100K_with_feature_size.csv", row.names=F, quote=F)

### 500000
All_genes_windowed_500K <- All_genes %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(All_genes_windowed_500K)
All_genes_windowed_500K<-All_genes_windowed_500K[,c(1:2,6:8,11)]
All_genes_windowed_500K<-All_genes_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)%>%mutate(window_start=(position_window*500000+1))
write.csv(All_genes_windowed_500K, "6_2_All_genes_windowed_500K_with_feature_size.csv", row.names=F, quote=F)

### 1000000
All_genes_windowed_1000K <- All_genes %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%1000000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(All_genes_windowed_1000K)
All_genes_windowed_1000K<-All_genes_windowed_1000K[,c(1:2,6:8,11)]
All_genes_windowed_1000K<-All_genes_windowed_1000K %>% mutate(proportion_of_window=featuresize_sum/1000000)%>%mutate(window_start=(position_window*1000000+1))
#plot(All_genes_windowed_1000K$window_start, All_genes_windowed_1000K$LG_start_mean)
write.csv(All_genes_windowed_1000K, "6_2_All_genes_windowed_1000K_with_feature_size.csv", row.names=F, quote=F)

#View(All_genes_windowed_500K)

getwd()

#write.csv(All_genes_windowed_500K, "All_genes_windowed_500K.csv")


BUSCO_genes_windowed_500K <- All_genes_BUSCO %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 

BUSCO_genes_windowed_500K <- All_genes_BUSCO %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(BUSCO_genes_windowed_500K)
BUSCO_genes_windowed_500K<-BUSCO_genes_windowed_500K[,c(1:2,10,20)]
BUSCO_genes_windowed_500K<-BUSCO_genes_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)


#write.csv(BUSCO_genes_windowed_500K, "BUSCO_genes_windowed_500K_with_proportion_including_repeat_overlap.csv")
#write.csv(BUSCO_genes_windowed_500K, "BUSCO_genes_windowed_500K.csv")



window_size<-50000
window_size<-100000
window_size<-500000
window_size<-1000000
window_size<-5000000


assign(paste("Repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("Repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("All_repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)



window_sizes<-c(50000, 100000, 500000, 1000000, 5000000)


for (window_size in window_sizes) {


assign(paste("MULE_Repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% filter(., transcript_status=="DNA/MuLE-MuDR") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("MULE_Repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("MULE_Repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)

}

for (window_size in window_sizes) {

assign(paste("LINE_Repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% filter(., transcript_status=="LINE/L1") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("LINE_Repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("LINE_Repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)

}


for (window_size in window_sizes) {

assign(paste("COPIA_Repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% filter(., transcript_status=="LTR/Copia") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("COPIA_Repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("COPIA_Repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)

}


for (window_size in window_sizes) {

assign(paste("Gypsy_Repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% filter(., transcript_status=="LTR/Gypsy") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("Gypsy_Repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("Gypsy_Repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)

}


for (window_size in window_sizes) {

assign(paste("Simple_repeats_windowed_",as.character(window_size),sep=""), All_repeats %>% filter(., transcript_status=="Simple_repeat") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>%  add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) %>% mutate(proportion_of_window=featuresize_sum/as.numeric(window_size)) %>% mutate (window_start=(position_window*(as.numeric(window_size)+1)))) 

current_repeats<-select(get ( paste("Simple_repeats_windowed_",as.character(window_size),sep="")), c( "LG", "position_window", "window_start", "n_mean", "proportion_of_window","featuresize_mean") )

write.csv(current_repeats, paste("Simple_repeats_windowed_",as.character(window_size), "_", Sys.Date(), ".csv", sep=""), row.names=F, quote=F)

}


##################################

MULE_windowed_500K <- All_repeats %>% filter(., transcript_status=="DNA/MuLE-MuDR") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 
write.csv(MULE_windowed_500K, "MULE_windowed_500K.csv")





LINE_windowed_500K <- All_repeats %>% filter(., transcript_status=="LINE/L1") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 
write.csv(LINE_windowed_500K, "LINE_windowed_500K.csv")

LINE_windowed_500K <- All_repeats  %>% filter(., transcript_status=="LINE/L1") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(LINE_windowed_500K)
LINE_windowed_500K<-LINE_windowed_500K[,c(1:2,10,20)]
LINE_windowed_500K<-LINE_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)
write.csv(LINE_windowed_500K, "LINE_windowed_500K_with_feature_size.csv")




Copia_windowed_500K <- All_repeats %>% filter(., transcript_status=="LTR/Copia") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 
write.csv(Copia_windowed_500K, "Copia_windowed_500K.csv")


Copia_windowed_500K <- All_repeats  %>% filter(., transcript_status=="LTR/Copia") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(Copia_windowed_500K)
Copia_windowed_500K<-Copia_windowed_500K[,c(1:2,10,20)]
Copia_windowed_500K<-Copia_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)
write.csv(Copia_windowed_500K, "Copia_windowed_500K_with_feature_size.csv")




Gypsy_windowed_500K <- All_repeats %>% filter(., transcript_status=="LTR/Gypsy") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 
write.csv(Gypsy_windowed_500K, "Gypsy_windowed_500K.csv")

Gypsy_windowed_500K <- All_repeats  %>% filter(., transcript_status=="LTR/Gypsy") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(Gypsy_windowed_500K)
Gypsy_windowed_500K<-Gypsy_windowed_500K[,c(1:2,10,20)]
Gypsy_windowed_500K<-Gypsy_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)
write.csv(Gypsy_windowed_500K, "Gypsy_windowed_500K_with_feature_size.csv")



Simple_repeat_windowed_500K <- All_repeats %>% filter(., transcript_status=="Simple_repeat") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% tally() 
write.csv(Simple_repeat_windowed_500K, "Simple_repeat_windowed_500K.csv")


Simple_repeat_windowed_500K <- All_repeats  %>% filter(., transcript_status=="Simple_repeat") %>% group_by(LG)  %>% mutate(position_window=LG_lower_start%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
colnames(Simple_repeat_windowed_500K)
Simple_repeat_windowed_500K<-Simple_repeat_windowed_500K[,c(1:2,10,20)]
Simple_repeat_windowed_500K<-Simple_repeat_windowed_500K %>% mutate(proportion_of_window=featuresize_sum/500000)
write.csv(Simple_repeat_windowed_500K, "Simple_repeat_windowed_500K_with_feature_size.csv")




plot(All_genes_windowed_100K$position_window,All_genes_windowed_100K$n, ylim=c(0,100))
plot(L.10$position_window,L.10$n, ylim=c(0,150))

plot(L.3$position_window,L.3$n, ylim=c(0,800))

L.10<-subset(All_genes_windowed_100K, All_genes_windowed_100K$LG=="L.10")
L.10<-subset(All_genes_windowed_500K, All_genes_windowed_500K$LG=="L.10")
L.10<-subset(BUSCO_genes_windowed_500K, BUSCO_genes_windowed_500K$LG=="L.10")
L.10<-subset(MULE_windowed_500K, MULE_windowed_500K$LG=="L.10")
L.10<-subset(LINE_windowed_500K, LINE_windowed_500K$LG=="L.10")
L.10<-subset(Copia_windowed_500K, Copia_windowed_500K$LG=="L.10")
L.10<-subset(Gypsy_windowed_500K, Gypsy_windowed_500K$LG=="L.10")
L.10<-subset(Simple_repeat_windowed_500K, Simple_repeat_windowed_500K$LG=="L.10")
L.3<-subset(Repeats_windowed_500K, Repeats_windowed_500K$LG=="L.3")

View(All_genes_windowed_100K)

unique(All_genes_windowed_10K$position_window)
View(All_genes_windowed_10K)

%>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T))) 
write.csv(Windowed_summary_lengths10K, "10KWindows.csv")
