library(qtl2)
library(qtl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
#sessionInfo()
#?cbind.scan1





#TX Import data from file ------------
rumexTXF2<- read_cross2 ("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/QTL/rumex_trimmed.yaml")
rumexTXF2<- read_cross2 ("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/QTL/rumex_trimmed.yaml")



Xcovar <- get_x_covar(iron)

rumexTXF2

### TX add pseudomarkers to map #####
map <- insert_pseudomarkers(rumexTXF2$gmap, step=1)
#pull.map(ipomoea$gmap)
pr <- calc_genoprob(rumexTXF2, map, error_prob=0.002)
#pr <- calc_genoprob(ipomoea, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)

### TX scan for QTL ###### 
out_binary <- scan1(pr, rumexTXF2$pheno[,1], model="binary", maxit=100000)
out_binary <- scan1(pr, rumexTXF2$pheno[,1], model="binary", maxit=1000, bintol = 1e-6, eta_max=19) #No idea if these are sensible values for bintol and eta_max - eta_max seems to have more of an influence. Converges up to eta_max = 19. 

#out <- scan1(pr, rumexTXF2$pheno[,1])

binarypeaks<-find_peaks(out_binary, map, threshold=3, drop = 1.5, expand2markers = T)

lodindex lodcolumn  chr      pos      lod    ci_lo   ci_hi
1        1    pheno1 L.10 74.38434 26.45222 55.48413 75.5209
#116369765 #CI LOW
#148948683 #CI HIGH AND PEAK

######### TX make plots ############

plot(out_binary, map)
plot(out_binary, map, chr = "L.10")
plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))



############ TX Permutation analysis
setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/QTL")
operm_bin <- scan1perm(pr, rumexTXF2$pheno[,1], model="binary",
                       n_perm=1000)
write.table(operm_bin, "permutation_analysis_12-17-2019.txt")
summary(operm_bin, alpha=c(0.05,0.01,0.001,0.0001))

> summary(operm_bin, alpha=c(0.05,0.01,0.001,0.0001))
LOD thresholds (1000 permutations)
pheno1
0.05    3.52
0.01    4.14
0.001   5.10
1e-04   6.18

write.csv(out_binary, "sex_scan.csv")
sex_scan<-read.csv("sex_scan.csv")
sex_scan<-sex_scan%>%separate(.,X, into = c ("scaffold", "scaffpos"))

## TX add LG positions to map #####

map_file_chromo<-read.csv("Post_chromonomer_final_edited_markers_positions_colocated_filled.csv")
map_file_chromo<-map_file_chromo %>% separate(., marker, into = c("scaffold", "scaffpos","LGpos"))

map_with_LODs<-left_join(map_file_chromo, sex_scan, by=c("scaffold","scaffpos"))

write.csv(map_with_LODs, "Sex_QTL_scan_output.csv")

LG_10<-subset(map_with_LODs, map_with_LODs$chr=="L.10")
plot(LG_10$LGpos, LG_10$pheno1)














#NC Import data from file ------------
rumexNCF2<- read_cross2 ("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/QTL/rumex_XXY_trimmed.yaml")

### NC add pseudomarkers to map #####
map <- insert_pseudomarkers(rumexNCF2$gmap, step=1)
#pull.map(ipomoea$gmap)
pr <- calc_genoprob(rumexNCF2, map, error_prob=0.002)
#pr <- calc_genoprob(ipomoea, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)

### NC scan for QTL ###### 
out_binary <- scan1(pr, rumexNCF2$pheno[,1], model="binary", maxit=1000, bintol = 1e-6, eta_max=19) #No idea if these are sensible values for bintol and eta_max - eta_max seems to have more of an influence. Converges up to eta_max = 19. 

#out <- scan1(pr, rumexTXF2$pheno[,1])

binarypeaks<-find_peaks(out_binary, map, threshold=3, drop = 1.5, expand2markers = T)
#Separated fusion
lodindex lodcolumn  chr pos      lod    ci_lo    ci_hi
1        1    pheno1 L.10  79 20.36489 65.58134 83.84658

#116369765 #CI LOW
#148948683 #CI HIGH AND PEAK

######### NC make plots ############

plot(out_binary, map)
plot(out_binary, map, chr = "L.10")
plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))



############ NC Permutation analysis
setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/QTL")
operm_bin <- scan1perm(pr, rumexNCF2$pheno[,1], model="binary",
                       n_perm=1000)
write.table(operm_bin, "NC_permutation_analysis_12-19-2019.txt")
summary(operm_bin, alpha=c(0.05,0.01,0.001,0.0001))

> summary(operm_bin, alpha=c(0.05,0.01,0.001,0.0001))
LOD thresholds (1000 permutations)
pheno1
0.05    3.53
0.01    4.45
0.001   5.68
1e-04   6.49

write.csv(out_binary, "sex_scan_NC.csv")
sex_scan<-read.csv("sex_scan_NC.csv")
sex_scan<-sex_scan%>%separate(.,X, into = c ("scaffold", "scaffpos"))
View(out_binary)
## NC add LG positions to map #####


map_file_chromo<-read.csv("NC_Post_chromonomer_final_edited_markers_positions_colocated_filled.csv")

map_file_chromo<-map_file_chromo %>% separate(., marker, into = c("scaffold", "scaffpos","LGpos"))

map_with_LODs<-left_join(map_file_chromo, sex_scan, by=c("scaffold","scaffpos"))

write.csv(map_with_LODs, "NC_Sex_QTL_scan_output.csv")

LG_10<-subset(map_with_LODs, map_with_LODs$chr=="L.10")
plot(LG_10$LGpos, LG_10$pheno1)
LG_5<-subset(map_with_LODs, map_with_LODs$chr=="L.5")
plot(LG_5$LGpos, LG_5$pheno1)






