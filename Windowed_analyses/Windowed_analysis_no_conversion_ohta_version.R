library(dplyr)
library(tidyr)
#library(ggplot2)
#library(vcfR)
#library(zoo)
#library(pkgbuild)
#library(slide)
#library(remotes)
#remotes::install_github("DavisVaughan/slide")
#install.packages("broom")
#library(broom)
options(scipen=999)


#setwd("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses")
#setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses")

setwd("/ohta/joanna.rifkin/HiCSNPData/Windowed_analyses")
getwd()

#TX_LGs<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/TX_LGs.csv")
#TX_LGs<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/TX_LGs.csv")




#*************** SNP collecting ***************#
################ Collect sex-linked loci and heterozygosity from Josh's crosses #################
##### Import and transform data #####


#cross012<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/MergedCrossPopAnalysis/Complete012AllCrossPopMerged_Q20_GQ20_All_NoInv.012", header = T, na.strings = "-1", stringsAsFactors = F)


#cross012<-read.table("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/Sex_linked_SNPs/Complete012JoshCrossFiltered_Q20_GQ20_All_NoInv.012", header = T, na.strings = "-1", stringsAsFactors = F)
#cross012<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/Sex_linked_SNPs/Complete012JoshCrossFiltered_Q20_GQ20_All_NoInv.012", header = T, na.strings = "-1", stringsAsFactors = F)
cross012<-read.table("/ohta/joanna.rifkin/HiCSNPData/JoshCross/Complete012JoshCrossFiltered_Q20_GQ20_All_NoInv.012")
head(cross012)

##View(cross012)
# Add column with scaffold as number
str(cross012)
cross012 <- separate(cross012, CHROM, into=c("ScnbKXS", "Scaffold", "HRSCAF"), sep = "[_//;]", remove = FALSE)
cross012$Scaffold<-as.numeric(cross012$Scaffold)
colnames(cross012)[5]<-"Position"


cross012$Position<-as.numeric(cross012$Position)



# Adding male and female allele frequency columns. Removing sample1_Index.4.TXF1_F4_, who turned out in another analysis to have a different father.
cross012<-mutate(cross012,JXNCFAlleleFreq=((sample8_Index.5.NCMOM_ +
                                              sample2_Index_13.NCF1_F1_ +
                                              sample2_Index_14.NCF1_F2_ +
                                              sample2_Index_15.NCF1_F3_ +
                                              sample2_Index_16.NCF1_F4_ +
                                              sample2_Index_18.NCF1_F5_ +
                                              sample2_Index_19.NCF1_F6_))/14)

cross012<-mutate(cross012,JXNCMAlleleFreq=((sample7_Index.11.NCDAD_ +
                                              sample2_Index_20.NCF1_M1_ +
                                              sample2_Index_21.NCF1_M2_ +
                                              sample2_Index_22.NCF1_M3_ +
                                              sample2_Index_23.NCF1_M4_ +
                                              sample2_Index_25.NCF1_M5_ +
                                              sample2_Index_27.NCF1_M6_ ))/14)

cross012<-mutate(cross012,JXTXMAlleleFreq=((sampleTXDAD_ +
                                              sample1_Index.7.TXF1_M1_ +
                                              sample1_Index.8.TXF1_M2_ +
                                              sample1_Index.9.TXF1_M3_ +
                                              sample1_Index.10.TXF1_M4_ +
                                              sample1_Index.11.TXF1_M5_ +
                                              sample1_Index.12.TXF1_M6_ ))/14)

cross012<-mutate(cross012,JXTXFAlleleFreq=((sampleTXMOM_ +
                                              sample1_Index.1.TXF1_F1_ +
                                              sample1_Index.2.TXF1_F2_ +
                                              sample1_Index.3.TXF1_F3_ +
                                              #sample1_Index.4.TXF1_F4_ +
                                              sample1_Index.5.TXF1_F5_ +
                                              sample1_Index.6.TXF1_F6_ ))/12)



# Adding counts of homozygous and heterozygous offspring. A = reference, B = nonreference
homozygousA<-c(0)
homozygousB<-c(2)
heterozygous<-c(1)
colnames(cross012)

cross012$JXTXSonsA<-rowSums((cross012[c(6:8,14:16)]==0)==T)
cross012$JXTXSonsB<-rowSums((cross012[c(6:8,14:16)]==2)==T)
cross012$JXTXSonsAB<-rowSums((cross012[c(6:8,14:16)]==1)==T)

cross012$JXTXDaughtersA<-rowSums((cross012[c(9:13)]==0)==T)
cross012$JXTXDaughtersB<-rowSums((cross012[c(9:13)]==2)==T)
cross012$JXTXDaughtersAB<-rowSums((cross012[c(9:13)]==1)==T)

cross012$JXNCSonsA<-rowSums((cross012[c(23:28)]==0)==T)
cross012$JXNCSonsB<-rowSums((cross012[c(23:28)]==2)==T)
cross012$JXNCSonsAB<-rowSums((cross012[c(23:28)]==1)==T)

cross012$JXNCDaughtersA<-rowSums((cross012[c(17:22)]==0)==T)
cross012$JXNCDaughtersB<-rowSums((cross012[c(17:22)]==2)==T)
cross012$JXNCDaughtersAB<-rowSums((cross012[c(17:22)]==1)==T)


cross012$JXTXFemaleHeterozygosity=cross012$JXTXDaughtersAB/(cross012$JXTXDaughtersAB + cross012$JXTXDaughtersA + cross012$JXTXDaughtersB)
cross012$JXTXMaleHeterozygosity=cross012$JXTXSonsAB/(cross012$JXTXSonsAB + cross012$JXTXSonsB + cross012$JXTXSonsA)

hist(cross012$JXTXFemaleHeterozygosity)
hist(cross012$JXTXMaleHeterozygosity)

colnames(cross012)
  ############# TX ##################


    ########## the presence of segregating Y-linked variants #########

homozygous<-c(0,2)
#homozygous<-c(0)
#homozygous<-c(2)
homozygousA<-c(0)
homozygousB<-c(2)
heterozygous<-c(1)
colnames(cross012)


#the presence of a segregating Y-linked variant, where fathers and sons were heterozygous but mothers and daughters were homozygous
#currently not broken down by homozygous ref / alt
cross012<-   cross012 %>%
  mutate(
    JXTX_YLinked = case_when(
      (
        sampleTXDAD_ %in% heterozygous & 
          sampleTXMOM_ %in% homozygous &
          sample1_Index.7.TXF1_M1_ %in% heterozygous & 
          sample1_Index.8.TXF1_M2_ %in% heterozygous &
          sample1_Index.9.TXF1_M3_ %in% heterozygous &
          sample1_Index.10.TXF1_M4_ %in% heterozygous &
          sample1_Index.11.TXF1_M5_ %in% heterozygous &
          sample1_Index.12.TXF1_M6_ %in% heterozygous &
          sample1_Index.1.TXF1_F1_ %in% homozygous &
          sample1_Index.2.TXF1_F2_ %in% homozygous &
          sample1_Index.3.TXF1_F3_ %in% homozygous &
          # sample1_Index.4.TXF1_F4_ %in% homozygous &
          sample1_Index.5.TXF1_F5_ %in% homozygous &
          sample1_Index.6.TXF1_F6_ %in% homozygous 
      ) ~ 1,
      TRUE ~ 0))

    ########## the presence of segregating X-linked variants #########
#the presence of a segregating X-linked variant, where fathers and daughters were heterozygous but mothers and sons were homozygous.


#homozygous<-c(0,2)#
#homozygous<-c(0)
#homozygous<-c(2)
#heterozygous<-c(1)


cross012<-   cross012 %>%
  mutate(
    JXTX_XLinked = case_when(
      (
        sampleTXDAD_ %in% heterozygous & 
          sampleTXMOM_ %in% homozygous &
          sample1_Index.7.TXF1_M1_ %in% homozygous & 
          sample1_Index.8.TXF1_M2_ %in% homozygous &
          sample1_Index.9.TXF1_M3_ %in% homozygous &
          sample1_Index.10.TXF1_M4_ %in% homozygous &
          sample1_Index.11.TXF1_M5_ %in% homozygous &
          sample1_Index.12.TXF1_M6_ %in% homozygous &
          sample1_Index.1.TXF1_F1_ %in% heterozygous &
          sample1_Index.2.TXF1_F2_ %in% heterozygous &
          sample1_Index.3.TXF1_F3_ %in% heterozygous &
          #sample1_Index.4.TXF1_F4_ %in% heterozygous &
          sample1_Index.5.TXF1_F5_ %in% heterozygous &
          sample1_Index.6.TXF1_F6_ %in% heterozygous
      ) ~ 1,
      TRUE ~ 0))

    ########## Hemizygous sites #########
homozygous<-c(0,2)
homozygousA<-c(0)
homozygousB<-c(2)
heterozygous<-c(1)


#a1) Maternal genotype AA, paternal genotype called BB. All daughters AB, all sons called AA
#a2) Maternal genotype BB, paternal genotype called AA. All daughters AB, all sons called BB
#b)Maternal genotype AB, paternal genotype called AA (or BB), some daughters AB, no sons heterozygous, the set of sons showing BOTH AA and BB calls 

cross012<-   cross012 %>%
  mutate(
    JXTX_Hemizygous = case_when(
      (
        sampleTXDAD_ %in% homozygousB & 
          sampleTXMOM_ %in% homozygousA &
          sample1_Index.7.TXF1_M1_ %in% homozygousA & 
          sample1_Index.8.TXF1_M2_ %in% homozygousA &
          sample1_Index.9.TXF1_M3_ %in% homozygousA &
          sample1_Index.10.TXF1_M4_ %in% homozygousA &
          sample1_Index.11.TXF1_M5_ %in% homozygousA &
          sample1_Index.12.TXF1_M6_ %in% homozygousA &
          sample1_Index.1.TXF1_F1_ %in% heterozygous &
          sample1_Index.2.TXF1_F2_ %in% heterozygous &
          sample1_Index.3.TXF1_F3_ %in% heterozygous &
          #sample1_Index.4.TXF1_F4_ %in% heterozygous &
          sample1_Index.5.TXF1_F5_ %in% heterozygous &
          sample1_Index.6.TXF1_F6_ %in% heterozygous
      ) ~ 1,
      (
        sampleTXDAD_ %in% homozygousA & 
          sampleTXMOM_ %in% homozygousB &
          sample1_Index.7.TXF1_M1_ %in% homozygousB & 
          sample1_Index.8.TXF1_M2_ %in% homozygousB &
          sample1_Index.9.TXF1_M3_ %in% homozygousB &
          sample1_Index.10.TXF1_M4_ %in% homozygousB &
          sample1_Index.11.TXF1_M5_ %in% homozygousB &
          sample1_Index.12.TXF1_M6_ %in% homozygousB &
          sample1_Index.1.TXF1_F1_ %in% heterozygous &
          sample1_Index.2.TXF1_F2_ %in% heterozygous &
          sample1_Index.3.TXF1_F3_ %in% heterozygous &
          #sample1_Index.4.TXF1_F4_ %in% heterozygous &
          sample1_Index.5.TXF1_F5_ %in% heterozygous &
          sample1_Index.6.TXF1_F6_ %in% heterozygous      
      )~ 1,
      (
        sampleTXDAD_ %in% homozygousB 
        &
          sampleTXMOM_ %in% heterozygous 
        &
          !(sample1_Index.7.TXF1_M1_ %in% heterozygous) &
          !(sample1_Index.8.TXF1_M2_ %in% heterozygous) &
          !(sample1_Index.9.TXF1_M3_ %in% heterozygous) &
          !(sample1_Index.10.TXF1_M4_ %in% heterozygous) &
          !(sample1_Index.11.TXF1_M5_ %in% heterozygous) &
          !(sample1_Index.12.TXF1_M6_ %in% heterozygous)
        &
          
          ((sample1_Index.7.TXF1_M1_ %in% homozygousB) |
             (sample1_Index.8.TXF1_M2_ %in% homozygousB) |
             (sample1_Index.9.TXF1_M3_ %in% homozygousB) |
             (sample1_Index.10.TXF1_M4_ %in% homozygousB) |
             (sample1_Index.11.TXF1_M5_ %in% homozygousB) |
             (sample1_Index.12.TXF1_M6_ %in% homozygousB))
        &
          ((sample1_Index.7.TXF1_M1_ %in% homozygousA) |
             (sample1_Index.8.TXF1_M2_ %in% homozygousA) |
             (sample1_Index.9.TXF1_M3_ %in% homozygousA) |
             (sample1_Index.10.TXF1_M4_ %in% homozygousA) |
             (sample1_Index.11.TXF1_M5_ %in% homozygousA) |
             (sample1_Index.12.TXF1_M6_ %in% homozygousA))
        &
          ((sample1_Index.1.TXF1_F1_ %in% heterozygous) |
             (sample1_Index.2.TXF1_F2_ %in% heterozygous )|
             (sample1_Index.3.TXF1_F3_ %in% heterozygous) |
             #sample1_Index.4.TXF1_F4_ %in% heterozygous &
             (sample1_Index.5.TXF1_F5_ %in% heterozygous) |
             (sample1_Index.6.TXF1_F6_ %in% heterozygous))
        
      )~ 1,
      TRUE ~ 0))

    ##### Autosomal sites ######
#To screen for autosomal genes, we identified SNPs where the mother was homozygous, father heterozygous, and at least two sons and at least two daughters were heterozygous (i.e. both sons and daughters inherit a focal allele from the father).

homozygous<-c(0,2)
#homozygous<-c(0,2)
#homozygousA<-c(0)
#homozygousB<-c(2)
heterozygous<-c(1)

colnames(cross012)
#cross012<-   cross012 %>%
#  mutate(
#    TX_Hemizygous = case_when(

cross012<-cross012 %>%
  mutate(
    JXTX_Autosomal = case_when(
      (                       sampleTXDAD_ %in% heterozygous 
                              &
                                sampleTXMOM_ %in% homozygous 
                              &
                                JXTXSonsAB > 1
                              & 
                                JXTXDaughtersAB > 1)~ 1,
      TRUE ~ 0))


##View(cross012)






############# NC ##################

colnames(cross012)



homozygous<-c(0,2)
#homozygous<-c(0)
#homozygous<-c(2)
homozygousA<-c(0)
homozygousB<-c(2)
heterozygous<-c(1)
colnames(cross012)

  ########## the presence of segregating Y-linked variants #########

#the presence of a segregating Y-linked variant, where fathers and sons were heterozygous but mothers and daughters were homozygous
#currently not broken down by homozygous ref / alt
cross012<-   cross012 %>%
  mutate(
    JXNC_YLinked = case_when(
      (
        sample7_Index.11.NCDAD_ %in% heterozygous & 
          sample8_Index.5.NCMOM_ %in% homozygous &
          sample2_Index_13.NCF1_F1_ %in% homozygous &
          sample2_Index_14.NCF1_F2_ %in% homozygous &
          sample2_Index_15.NCF1_F3_ %in% homozygous &
          sample2_Index_16.NCF1_F4_ %in% homozygous &
          sample2_Index_18.NCF1_F5_ %in% homozygous &
          sample2_Index_19.NCF1_F6_ %in% homozygous &
          sample2_Index_20.NCF1_M1_ %in% heterozygous & 
          sample2_Index_21.NCF1_M2_ %in% heterozygous &
          sample2_Index_22.NCF1_M3_ %in% heterozygous &
          sample2_Index_23.NCF1_M4_ %in% heterozygous &
          sample2_Index_25.NCF1_M5_ %in% heterozygous &
          sample2_Index_27.NCF1_M6_ %in% heterozygous 
      ) ~ 1,
      TRUE ~ 0))


  ########## the presence of segregating X-linked variants #########
#the presence of a segregating X-linked variant, where fathers and daughters were heterozygous but mothers and sons were homozygous.

cross012<-   cross012 %>%
  mutate(
    JXNC_XLinked = case_when(
      (
        sample7_Index.11.NCDAD_ %in% heterozygous & 
          sample8_Index.5.NCMOM_ %in% homozygous &
          sample2_Index_13.NCF1_F1_ %in% heterozygous &
          sample2_Index_14.NCF1_F2_ %in% heterozygous &
          sample2_Index_15.NCF1_F3_ %in% heterozygous &
          sample2_Index_16.NCF1_F4_ %in% heterozygous &
          sample2_Index_18.NCF1_F5_ %in% heterozygous &
          sample2_Index_19.NCF1_F6_ %in% heterozygous &
          sample2_Index_20.NCF1_M1_ %in% homozygous & 
          sample2_Index_21.NCF1_M2_ %in% homozygous &
          sample2_Index_22.NCF1_M3_ %in% homozygous &
          sample2_Index_23.NCF1_M4_ %in% homozygous &
          sample2_Index_25.NCF1_M5_ %in% homozygous &
          sample2_Index_27.NCF1_M6_ %in% homozygous 
      ) ~ 1,
      TRUE ~ 0))

homozygous<-c(0,2)
homozygousA<-c(0)
homozygousB<-c(2)
heterozygous<-c(1)

  ########## Hemizygous sites #########
#a1) Maternal genotype AA, paternal genotype called BB. All daughters AB, all sons called AA
#a2) Maternal genotype BB, paternal genotype called AA. All daughters AB, all sons called BB
#b)Maternal genotype AB, paternal genotype called AA (or BB), some daughters AB, no sons heterozygous, the set of sons showing BOTH AA and BB calls 

cross012<-   cross012 %>%
  mutate(
    JXNC_Hemizygous = case_when(
      (
        sample7_Index.11.NCDAD_ %in% homozygousB & 
          sample8_Index.5.NCMOM_ %in% homozygousA &
          
          sample2_Index_20.NCF1_M1_ %in% homozygousA & 
          sample2_Index_21.NCF1_M2_ %in% homozygousA &
          sample2_Index_22.NCF1_M3_ %in% homozygousA &
          sample2_Index_23.NCF1_M4_ %in% homozygousA &
          sample2_Index_25.NCF1_M5_ %in% homozygousA &
          sample2_Index_27.NCF1_M6_ %in% homozygousA &
          
          sample2_Index_13.NCF1_F1_ %in% heterozygous &
          sample2_Index_14.NCF1_F2_ %in% heterozygous &
          sample2_Index_15.NCF1_F3_ %in% heterozygous &
          sample2_Index_16.NCF1_F4_ %in% heterozygous &
          sample2_Index_18.NCF1_F5_ %in% heterozygous &
          sample2_Index_19.NCF1_F6_ %in% heterozygous
      ) ~ 1,
      (
        sample7_Index.11.NCDAD_ %in% homozygousA & 
          sample8_Index.5.NCMOM_ %in% homozygousB &
          
          sample2_Index_20.NCF1_M1_ %in% homozygousB & 
          sample2_Index_21.NCF1_M2_ %in% homozygousB &
          sample2_Index_22.NCF1_M3_ %in% homozygousB &
          sample2_Index_23.NCF1_M4_ %in% homozygousB &
          sample2_Index_25.NCF1_M5_ %in% homozygousB &
          sample2_Index_27.NCF1_M6_ %in% homozygousB &
          
          sample2_Index_13.NCF1_F1_ %in% heterozygous &
          sample2_Index_14.NCF1_F2_ %in% heterozygous &
          sample2_Index_15.NCF1_F3_ %in% heterozygous &
          sample2_Index_16.NCF1_F4_ %in% heterozygous &
          sample2_Index_18.NCF1_F5_ %in% heterozygous &
          sample2_Index_19.NCF1_F6_ %in% heterozygous     
      )~ 1,
      (
        sample7_Index.11.NCDAD_ %in% homozygousB 
        &
          sample8_Index.5.NCMOM_ %in% heterozygous 
        &
          !(sample2_Index_20.NCF1_M1_ %in% heterozygous) &
          !(sample2_Index_21.NCF1_M2_ %in% heterozygous) &
          !(sample2_Index_22.NCF1_M3_ %in% heterozygous) &
          !(sample2_Index_23.NCF1_M4_ %in% heterozygous) &
          !(sample2_Index_25.NCF1_M5_ %in% heterozygous) &
          !(sample2_Index_27.NCF1_M6_ %in% heterozygous)
        &
          
          ((sample2_Index_20.NCF1_M1_ %in% homozygousB) |
             (sample2_Index_21.NCF1_M2_ %in% homozygousB) |
             (sample2_Index_22.NCF1_M3_ %in% homozygousB) |
             (sample2_Index_23.NCF1_M4_ %in% homozygousB) |
             (sample2_Index_25.NCF1_M5_ %in% homozygousB) |
             (sample2_Index_27.NCF1_M6_ %in% homozygousB))
        &
          ((sample2_Index_20.NCF1_M1_ %in% homozygousA) |
             (sample2_Index_21.NCF1_M2_ %in% homozygousA) |
             (sample2_Index_22.NCF1_M3_ %in% homozygousA) |
             (sample2_Index_23.NCF1_M4_ %in% homozygousA) |
             (sample2_Index_25.NCF1_M5_ %in% homozygousA) |
             (sample2_Index_27.NCF1_M6_ %in% homozygousA))
        &
          ((sample2_Index_13.NCF1_F1_ %in% heterozygous) |
             (sample2_Index_14.NCF1_F2_ %in% heterozygous )|
             (sample2_Index_15.NCF1_F3_ %in% heterozygous) |
             (sample2_Index_16.NCF1_F4_ %in% heterozygous) |
             (sample2_Index_18.NCF1_F5_ %in% heterozygous) |
             (sample2_Index_19.NCF1_F6_ %in% heterozygous))
        
      )~ 1,
      TRUE ~ 0))

  ##### Autosomal sites ######

cross012<-cross012 %>%
  mutate(
    JXNC_Autosomal = case_when(
      (                        sample7_Index.11.NCDAD_ %in% heterozygous 
                               &
                                 sample8_Index.5.NCMOM_ %in% homozygous 
                               &
                                 JXNCSonsAB > 1
                               & 
                                 JXNCDaughtersAB > 1)~ 1,
      TRUE ~ 0))

######## Shared pedigree SNPs NC and TX ######

colnames(cross012)

cross012<-cross012 %>%
  mutate(
    JX_shared_X_linked = case_when(
      (                        JXTX_XLinked == 1
                               &
                                 JXNC_XLinked == 1)~ 1,
      TRUE ~ 0))


cross012<-cross012 %>%
  mutate(
    JX_shared_Y_linked = case_when(
      (                        JXTX_YLinked == 1
                               &
                                 JXNC_YLinked == 1)~ 1,
      TRUE ~ 0))

cross012<-cross012 %>%
  mutate(
    JX_shared_Hemizygous = case_when(
      (                        JXTX_Hemizygous == 1
                               &
                                 JXNC_Hemizygous == 1)~ 1,
      TRUE ~ 0))

cross012<-cross012 %>%
  mutate(
    JX_shared_Autosomal = case_when(
      (                        JXTX_Autosomal == 1
                               &
                                 JXNC_Autosomal == 1)~ 1,
      TRUE ~ 0))




##View(subset(cross012, cross012$JXTX_XLinked==1))

colnames(cross012)
pedigree_summary<-cross012[,c(1:5,77:106)]
rm(cross012)
# write.csv(pedigree_summary, "12-3-2019_pedigree_summary.csv") #Writing a backup file

# write.csv(cross012, "11-25-2019_cross_with_AF_and_Offspring_Alleles.csv") #Writing a backup file


rm(cross012)
################ F2 Collect sex-associated loci and heterozygosity from the F2 linkage mapping dataset #################





############# TX ##################
######-------------- Import data frame from 012 conversion pipeline ------
TTfile<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/PosIndvAddedCopyTX_GQ50_0.05_correct_ind_only.csvr", header=T, stringsAsFactors = F)
TTfile<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/PosIndvAddedCopyTX_GQ50_0.05_correct_ind_only.csvr", header=T, stringsAsFactors = F)
##View(TTfile)


TTfile <- separate(TTfile, Genotype, c("ScnbKXS","Scaffold", "HRSCAF", "HRSCAFnum","Position"), remove=FALSE)
TTfile$Scaffold<-as.numeric(TTfile$Scaffold)
TTfile$Position<-as.numeric(TTfile$Position)


############## TX allele frequencies ###########

TTfile$F2TXMalesA<-rowSums((TTfile[c(9,14,15,16,19,20,22,23,25,26,27,28,30,32,35,37,38,40,43,45,47,50,51,52,53,55,56,58,60,61,62,63,64,67,69,75,76,77,78,82,83,86,88,91,92,93,96,97,98,99,100,101)]=="A")==T)
TTfile$F2TXMalesB<-rowSums((TTfile[c(9,14,15,16,19,20,22,23,25,26,27,28,30,32,35,37,38,40,43,45,47,50,51,52,53,55,56,58,60,61,62,63,64,67,69,75,76,77,78,82,83,86,88,91,92,93,96,97,98,99,100,101)]=="B")==T)
TTfile$F2TXMalesH<-rowSums((TTfile[c(9,14,15,16,19,20,22,23,25,26,27,28,30,32,35,37,38,40,43,45,47,50,51,52,53,55,56,58,60,61,62,63,64,67,69,75,76,77,78,82,83,86,88,91,92,93,96,97,98,99,100,101)]=="H")==T)
TTfile$F2TXMaleHeterozygosity=as.numeric(TTfile$F2TXMalesH)/as.numeric(TTfile$F2TXMalesA+TTfile$F2TXMalesB+TTfile$F2TXMalesH)

TTfile$F2TXFemalesA<-rowSums((TTfile[c(8,10,11,12,13,17,18,21,24,29,31,33,34,36,39,41,42,44,46,48,49,54,57,59,65,66,68,70,71,72,73,74,79,80,81,84,85,87,89,90,94,95)]=="A")==T)
TTfile$F2TXFemalesB<-rowSums((TTfile[c(8,10,11,12,13,17,18,21,24,29,31,33,34,36,39,41,42,44,46,48,49,54,57,59,65,66,68,70,71,72,73,74,79,80,81,84,85,87,89,90,94,95)]=="B")==T)
TTfile$F2TXFemalesH<-rowSums((TTfile[c(8,10,11,12,13,17,18,21,24,29,31,33,34,36,39,41,42,44,46,48,49,54,57,59,65,66,68,70,71,72,73,74,79,80,81,84,85,87,89,90,94,95)]=="H")==T)
TTfile$F2TXFemaleHeterozygosity=as.numeric(TTfile$F2TXFemalesH)/as.numeric(TTfile$F2TXFemalesA+TTfile$F2TXFemalesB+TTfile$F2TXFemalesH)

hist(TTfile$F2TXFemaleHeterozygosity)
hist(TTfile$F2TXMaleHeterozygosity)

TTfile$F2TXAF_A=as.numeric(TTfile$F2TXFemalesH/2+TTfile$F2TXMalesH/2+TTfile$F2TXFemalesA+TTfile$F2TXMalesA)/as.numeric(TTfile$F2TXFemalesA+TTfile$F2TXFemalesB+TTfile$F2TXFemalesH + TTfile$F2TXMalesA+TTfile$F2TXMalesB+TTfile$F2TXMalesH)
hist(TTfile$F2TXAF_A)
TTfile$F2TXAF_B=as.numeric(TTfile$F2TXFemalesH/2+TTfile$F2TXMalesH/2+TTfile$F2TXFemalesB+TTfile$F2TXMalesB)/as.numeric(TTfile$F2TXFemalesA+TTfile$F2TXFemalesB+TTfile$F2TXFemalesH + TTfile$F2TXMalesA+TTfile$F2TXMalesB+TTfile$F2TXMalesH)
hist(TTfile$F2TXAF_B)


TTfile$F2TXMaleAF_A=as.numeric(TTfile$F2TXMalesH/2+TTfile$F2TXMalesA)/as.numeric(
  TTfile$F2TXMalesA+TTfile$F2TXMalesB+TTfile$F2TXMalesH)
hist(TTfile$F2TXMaleAF_A)
TTfile$F2TXMaleAF_B=as.numeric(TTfile$F2TXMalesH/2+TTfile$F2TXMalesB)/as.numeric(
  TTfile$F2TXMalesA+TTfile$F2TXMalesB+TTfile$F2TXMalesH)
hist(TTfile$F2TXMaleAF_B)

head(TTfile)

################ TX XY and XX/X- polymorphisms ##########

TTfile$XY=F2TXMaleHeterozygosity>0.9 & F2TXFemaleHeterozygosity<0.1



TTfile<-   TTfile %>%
  mutate(
    F2_XY_pattern = case_when(
      (
        F2TXMaleHeterozygosity > 0.9 & 
          F2TXFemaleHeterozygosity < 0.1) ~ 1,
      TRUE ~ 0))

TTfile<-   TTfile %>%
  mutate(
    F2_Hemi_pattern = case_when(
      (
        F2TXFemaleHeterozygosity > 0.9 & 
          F2TXMaleHeterozygosity < 0.1) ~ 1,
      TRUE ~ 0))

sum(TTfile$F2_XY_pattern)
sum(TTfile$F2_Hemi_pattern)

############# NC ##################
######-------------- Import data frame from 012 conversion pipeline 
setwd("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/NC_transcriptome")


NCfile<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/NC_transcriptome/PosIndvAddedCopyNC_GQ50_0.05_correct_ind_only.csvr", header=T, stringsAsFactors = F)

NCfile<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/NC_transcriptome/PosIndvAddedCopyNC_GQ50_0.05_correct_ind_only.csvr", header=T, stringsAsFactors = F)
##View(TTfile)


NCfile <- separate(NCfile, Genotype, c("ScnbKXS","Scaffold", "HRSCAF", "HRSCAFnum","Position"), remove=FALSE)
NCfile$Scaffold<-as.numeric(NCfile$Scaffold)
NCfile$Position<-as.numeric(NCfile$Position)


colnames(NCfile)


############## NC allele frequencies ###########
NCfile$F2NCMalesA<-rowSums((NCfile[c(8,10,16,17,21,22,25,26,28,30,33,34,45,46,47,50,51,53,54,55,57,61,68,73,75,76,77,78,79,80,83,84,87,88,89,95,97,100)]=="A")==T)
NCfile$F2NCMalesB<-rowSums((NCfile[c(8,10,16,17,21,22,25,26,28,30,33,34,45,46,47,50,51,53,54,55,57,61,68,73,75,76,77,78,79,80,83,84,87,88,89,95,97,100)]=="B")==T)
NCfile$F2NCMalesH<-rowSums((NCfile[c(8,10,16,17,21,22,25,26,28,30,33,34,45,46,47,50,51,53,54,55,57,61,68,73,75,76,77,78,79,80,83,84,87,88,89,95,97,100)]=="H")==T)

NCfile$F2NCMaleHeterozygosity=as.numeric(NCfile$F2NCMalesH)/as.numeric(NCfile$F2NCMalesA+NCfile$F2NCMalesB+NCfile$F2NCMalesH)

NCfile$F2NCFemalesA<-rowSums((NCfile[c(4,9,11,12,13,14,15,18,19,20,23,24,27,29,31,32,35,36,37,38,39,40,41,42,43,44,48,49,52,56,58,59,60,62,63,64,65,66,67,69,70,71,72,74,81,82,85,86,90,91,92,93,94,96,98,99)]=="A")==T)
NCfile$F2NCFemalesB<-rowSums((NCfile[c(4,9,11,12,13,14,15,18,19,20,23,24,27,29,31,32,35,36,37,38,39,40,41,42,43,44,48,49,52,56,58,59,60,62,63,64,65,66,67,69,70,71,72,74,81,82,85,86,90,91,92,93,94,96,98,99)]=="B")==T)
NCfile$F2NCFemalesH<-rowSums((NCfile[c(4,9,11,12,13,14,15,18,19,20,23,24,27,29,31,32,35,36,37,38,39,40,41,42,43,44,48,49,52,56,58,59,60,62,63,64,65,66,67,69,70,71,72,74,81,82,85,86,90,91,92,93,94,96,98,99)]=="H")==T)
NCfile$F2NCFemaleHeterozygosity=as.numeric(NCfile$F2NCFemalesH)/as.numeric(NCfile$F2NCFemalesA+NCfile$F2NCFemalesB+NCfile$F2NCFemalesH)




NCfile$F2NCAF_A=as.numeric(NCfile$F2NCFemalesH/2+NCfile$F2NCMalesH/2+NCfile$F2NCFemalesA+NCfile$F2NCMalesA)/as.numeric(NCfile$F2NCFemalesA+NCfile$F2NCFemalesB+NCfile$F2NCFemalesH + NCfile$F2NCMalesA+NCfile$F2NCMalesB+NCfile$F2NCMalesH)
hist(NCfile$F2NCAF_A)
NCfile$F2NCAF_B=as.numeric(NCfile$F2NCFemalesH/2+NCfile$F2NCMalesH/2+NCfile$F2NCFemalesB+NCfile$F2NCMalesB)/as.numeric(NCfile$F2NCFemalesA+NCfile$F2NCFemalesB+NCfile$F2NCFemalesH + NCfile$F2NCMalesA+NCfile$F2NCMalesB+NCfile$F2NCMalesH)
hist(NCfile$F2NCAF_B)


NCfile$F2NCMaleAF_A=as.numeric(NCfile$F2NCMalesH/2+NCfile$F2NCMalesA)/as.numeric(
  NCfile$F2NCMalesA+NCfile$F2NCMalesB+NCfile$F2NCMalesH)
hist(NCfile$F2NCMaleAF_A)
NCfile$F2NCMaleAF_B=as.numeric(NCfile$F2NCMalesH/2+NCfile$F2NCMalesB)/as.numeric(
  NCfile$F2NCMalesA+NCfile$F2NCMalesB+NCfile$F2NCMalesH)
hist(NCfile$F2NCMaleAF_B)

NCfile<-   NCfile %>%
  mutate(
    F2_XY_pattern = case_when(
      (
        F2NCMaleHeterozygosity > 0.9 & 
          F2NCFemaleHeterozygosity < 0.1) ~ 1,
      TRUE ~ 0))

NCfile<-   NCfile %>%
  mutate(
    F2_Hemi_pattern = case_when(
      (
        F2NCFemaleHeterozygosity > 0.9 & 
          F2NCMaleHeterozygosity < 0.1) ~ 1,
      TRUE ~ 0))

sum(NCfile$F2_XY_pattern)
sum(NCfile$F2_Hemi_pattern)

########### store F2 data #############

head(NCfile1)
colnames(NCfile)
colnames(TTfile)

NCfile1<-NCfile[,c(1,3, 6:7, 101:114)]
NCfile1<-NCfile1[-1,]
TTfile1<-TTfile[,c(1,3, 6:7, 102:115)]
TTfile1<-TTfile1[-1,]


F2_summary<-full_join(NCfile1, TTfile1, by=c("Genotype", "Scaffold", "Position"))

write.csv(F2_summary, "12-3-2019_F2_summary.csv")

rm(NCfile)
rm(TTfile)

################ TX (Josh and Felix's data) ##################
################ NC (Josh and Felix's data) ##################

################ Collect fixed differences, heterozygosity and shared polymorphisms between NC and TX pop data ##################

cross012<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/Completefinal012JoshFelixPopFiltered_Q20_GQ20_All_NoInv.012", header = T, na.strings = "-1", stringsAsFactors = F)
cross012<-cross012[-1,] #Remove positio row

colnames(cross012)
head(cross012)
cross012 <- separate(cross012, CHROM, into=c("ScnbKXS", "Scaffold", "HRSCAF"), sep = "[_//;]", remove = FALSE)
colnames(all_missing)



################ Allele frequencies #################



cross012<-mutate(cross012,JPRefNCAlleleFreq=((sample20.NM1_ +
                                           sample21.NM2_ +
                                           sample22.NM3_ +
                                           sample23.NM4_ +
                                           sample25.NM5_ +
                                           sample27.NM6_ +
                                           sample13.NF1_ +
                                           sample14.NF2_ +
                                           sample15.NF3_ +
                                           sample16.NF4_ +
                                           sample18.NF5_ +
                                           sample19.NF6_))/24)


cross012<-mutate(cross012,JPRefNCMAlleleFreq=((sample20.NM1_ +
                                            sample21.NM2_ +
                                            sample22.NM3_ +
                                            sample23.NM4_ +
                                            sample25.NM5_ +
                                            sample27.NM6_))/12)

cross012<-mutate(cross012,JPRefNCFAlleleFreq=((sample13.NF1_ +
                                            sample14.NF2_ +
                                            sample15.NF3_ +
                                            sample16.NF4_ +
                                            sample18.NF5_ +
                                            sample19.NF6_))/12) 


cross012<-mutate(cross012,JPRefTXAlleleFreq=((sample7.TM1_  +
                                           sample8.TM2_ +
                                           sample9.TM3_ +
                                           sample10.TM4_ +
                                           sample11.TM5_ +
                                           sample12.TM6_ +
                                           sample1.TF1_ +
                                           sample2.TF2_ +
                                           sample3.TF3_ +
                                           sample4.TF4_ +
                                           sample5.TF5_ +
                                           sample6.TF6_))/24)


cross012<-mutate(cross012,JPRefTXMAlleleFreq=((sample7.TM1_ +
                                            sample8.TM2_ +
                                            sample9.TM3_ +
                                            sample10.TM4_ +
                                            sample11.TM5_ +
                                            sample12.TM6_))/12) 

cross012<-mutate(cross012,JPRefTXFAlleleFreq=((sample1.TF1_ +
                                            sample2.TF2_ +
                                            sample3.TF3_ +
                                            sample4.TF4_ +
                                            sample5.TF5_ +
                                            sample6.TF6_))/12)





cross012<-mutate(cross012,JPFPRefNCAlleleFreq=
                   ((sample20.NM1_ +
                     sample21.NM2_ +
                     sample22.NM3_ +
                     sample23.NM4_ +
                     sample25.NM5_ +
                       sample27.NM6_ +
                       sample13.NF1_ +
                       sample14.NF2_ +
                       sample15.NF3_ +
                       sample16.NF4_ +
                       sample18.NF5_ +
                       sample19.NF6_ +
                       sampleALBRE11 +
                       sampleALBRU8 +
                       sampleFLHAM4 +
                       sampleFLJAS13 +
                       sampleGABEL9 +
                       sampleGAGLA21 +
                       sampleGASTA3 +
                       sampleNCBAT4 +
                       sampleNCELI16 +
                       sampleNCHIC16 +
                       sampleNCKIN2 +
                       sampleNCROS7 +
                       sampleSCBRA4 +
                       sampleSCMAR24 +
                       sampleSCPRO27
                     ))/54)


cross012<-mutate(cross012,JPFPRefNCMAlleleFreq=
                   ((sample20.NM1_ +
                       sample21.NM2_ +
                       sample22.NM3_ +
                       sample23.NM4_ +
                       sample25.NM5_ +
                       sample27.NM6_ +
                       sampleALBRE11 +
                       sampleALBRU8 +
                       sampleFLHAM4 +
                       sampleFLJAS13 +
                       sampleGABEL9 +
                       sampleGAGLA21 +
                       sampleGASTA3 +
                       sampleNCBAT4 +
                       sampleNCELI16 +
                       sampleNCHIC16 +
                       sampleNCKIN2 +
                       sampleNCROS7 +
                       sampleSCBRA4 +
                       sampleSCMAR24 +
                       sampleSCPRO27))/42)

cross012<-mutate(cross012,JPFPRefTXAlleleFreq=
                   ((sample7.TM1_  +
                       sample8.TM2_ +
                       sample9.TM3_ +
                       sample10.TM4_ +
                       sample11.TM5_ +
                       sample12.TM6_ +
                       sample1.TF1_ +
                       sample2.TF2_ +
                       sample3.TF3_ +
                       sample4.TF4_ +
                       sample5.TF5_ +
                       sample6.TF6_+
                       sampleLABEN5 +
                       sampleOKBAC15 +
                       sampleOKRAT17 +
                       sampleTXATH5 +
                       sampleTXLIV14 +
                       sampleTXMTP16 +
                       sampleTXOAK6 +
                       sampleTXROS24))/40)


cross012<-mutate(cross012,JPFPRefTXMAlleleFreq=
                   ((sample7.TM1_ +
                       sample8.TM2_ +
                       sample9.TM3_ +
                       sample10.TM4_ +
                       sample11.TM5_ +
                       sample12.TM6_+
                       sampleLABEN5 +
                       sampleOKBAC15 +
                       sampleOKRAT17 +
                       sampleTXATH5 +
                       sampleTXLIV14 +
                       sampleTXMTP16 +
                       sampleTXOAK6 +
                       sampleTXROS24))/28) 



cross012<-mutate(cross012,FPRefNCAlleleFreq=
                   ((  sampleALBRE11 +
                         sampleALBRU8 +
                         sampleFLHAM4 +
                         sampleFLJAS13 +
                         sampleGABEL9 +
                         sampleGAGLA21 +
                         sampleGASTA3 +
                         sampleNCBAT4 +
                         sampleNCELI16 +
                         sampleNCHIC16 +
                         sampleNCKIN2 +
                         sampleNCROS7 +
                         sampleSCBRA4 +
                         sampleSCMAR24 +
                         sampleSCPRO27
                   ))/32)


cross012<-mutate(cross012,FPRefTXAlleleFreq=
                   ((sampleLABEN5 +
                       sampleOKBAC15 +
                       sampleOKRAT17 +
                       sampleTXATH5 +
                       sampleTXLIV14 +
                       sampleTXMTP16 +
                       sampleTXOAK6 +
                       sampleTXROS24))/16)


cross012$JP_NC_TX_Fixed<-as.numeric((cross012$JPRefTXAlleleFreq==0 & cross012$JPRefNCAlleleFreq==1) | (cross012$JPRefTXAlleleFreq==1 & cross012$JPRefNCAlleleFreq==0))

cross012$JP_NC_TX_Shared<-as.numeric((cross012$JPRefTXAlleleFreq>0 & cross012$JPRefTXAlleleFreq<1) & (cross012$JPRefNCAlleleFreq<1 & cross012$JPRefNCAlleleFreq>0))

cross012$JP_PolyNCOnly<-as.numeric((cross012$JPRefTXAlleleFreq==0 | cross012$JPRefTXAlleleFreq==1) & (cross012$JPRefNCAlleleFreq<1 & cross012$JPRefNCAlleleFreq>0)) 

cross012$JP_PolyTXOnly<-as.numeric((cross012$JPRefNCAlleleFreq==0 | cross012$JPRefNCAlleleFreq==1) & (cross012$JPRefTXAlleleFreq<1 & cross012$JPRefTXAlleleFreq>0)) 



cross012$FP_NC_TX_Fixed<-as.numeric((cross012$FPRefTXAlleleFreq==0 & cross012$FPRefNCAlleleFreq==1) | (cross012$FPRefTXAlleleFreq==1 & cross012$FPRefNCAlleleFreq==0))

cross012$FP_NC_TX_Shared<-as.numeric((cross012$FPRefTXAlleleFreq>0 & cross012$FPRefTXAlleleFreq<1) & (cross012$FPRefNCAlleleFreq<1 & cross012$FPRefNCAlleleFreq>0))

cross012$FP_PolyNCOnly<-as.numeric((cross012$FPRefTXAlleleFreq==0 | cross012$FPRefTXAlleleFreq==1) & (cross012$FPRefNCAlleleFreq<1 & cross012$FPRefNCAlleleFreq>0)) 

cross012$FP_PolyTXOnly<-as.numeric((cross012$FPRefNCAlleleFreq==0 | cross012$FPRefNCAlleleFreq==1) & (cross012$FPRefTXAlleleFreq<1 & cross012$FPRefTXAlleleFreq>0)) 




cross012$JPFP_NC_TX_Fixed<-as.numeric((cross012$JPFPRefTXAlleleFreq==0 & cross012$JPFPRefNCAlleleFreq==1) | (cross012$JPFPRefTXAlleleFreq==1 & cross012$JPFPRefNCAlleleFreq==0))

cross012$JPFP_NC_TX_Shared<-as.numeric((cross012$JPFPRefTXAlleleFreq>0 & cross012$JPFPRefTXAlleleFreq<1) & (cross012$JPFPRefNCAlleleFreq<1 & cross012$JPFPRefNCAlleleFreq>0))

cross012$JPFP_PolyNCOnly<-as.numeric((cross012$JPFPRefTXAlleleFreq==0 | cross012$JPFPRefTXAlleleFreq==1) & (cross012$JPFPRefNCAlleleFreq<1 & cross012$JPFPRefNCAlleleFreq>0)) 

cross012$JPFP_PolyTXOnly<-as.numeric((cross012$JPFPRefNCAlleleFreq==0 | cross012$JPFPRefNCAlleleFreq==1) & (cross012$JPFPRefTXAlleleFreq<1 & cross012$JPFPRefTXAlleleFreq>0)) 





colnames(cross012)
cross012$JP_TX_A<-rowSums((cross012[c(33:35,42,49:56)]==0)==T)
cross012$JP_TX_B<-rowSums((cross012[c(33:35,42,49:56)]==2)==T)
cross012$JP_TX_AB<-rowSums((cross012[c(33:35,42,49:56)]==1)==T)

cross012$JP_TXM_A<-rowSums((cross012[c(33:35,54:56)]==0)==T)
cross012$JP_TXM_B<-rowSums((cross012[c(33:35,54:56)]==2)==T)
cross012$JP_TXM_AB<-rowSums((cross012[c(33:35,54:56)]==1)==T)

cross012$JP_TXF_A<-rowSums((cross012[c(42,49:53)]==0)==T)
cross012$JP_TXF_B<-rowSums((cross012[c(42,49:53)]==2)==T)
cross012$JP_TXF_AB<-rowSums((cross012[c(42,49:53)]==1)==T)

cross012$FP_TX_A<-rowSums((cross012[c(33:35,42,49:56)]==0)==T)
cross012$FP_TX_B<-rowSums((cross012[c(33:35,42,49:56)]==2)==T)
cross012$FP_TX_AB<-rowSums((cross012[c(33:35,42,49:56)]==1)==T)

cross012$FP_TXM_A<-rowSums((cross012[c(33:35,54:56)]==0)==T)
cross012$FP_TXM_B<-rowSums((cross012[c(33:35,54:56)]==2)==T)
cross012$FP_TXM_AB<-rowSums((cross012[c(33:35,54:56)]==1)==T)


colnames(cross012)
cross012$FP_NC_A<-rowSums((cross012[c(36:48)]==0)==T)
cross012$FP_NC_B<-rowSums((cross012[c(36:48)]==2)==T)
cross012$FP_NC_AB<-rowSums((cross012[c(36:48)]==1)==T)

cross012$FP_NCM_A<-rowSums((cross012[c(43:48)]==0)==T)
cross012$FP_NCM_B<-rowSums((cross012[c(43:48)]==2)==T)
cross012$FP_NCM_AB<-rowSums((cross012[c(43:48)]==1)==T)


cross012$JPFP_TX_A<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,42,49:56)]==0)==T)
cross012$JPFP_TX_B<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,42,49:56)]==2)==T)
cross012$JPFP_TX_AB<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,42,49:56)]==1)==T)

cross012$JPFP_TXM_A<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,54:56)]==0)==T)
cross012$JPFP_TXM_B<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,54:56)]==2)==T)
cross012$JPFP_TXM_AB<-rowSums((cross012[c(10, 16:17, 21:25, 33:35,54:56)]==1)==T)





colnames(cross012)
cross012$JP_NC_A<-rowSums((cross012[c(36:48)]==0)==T)
cross012$JP_NC_B<-rowSums((cross012[c(36:48)]==2)==T)
cross012$JP_NC_AB<-rowSums((cross012[c(36:48)]==1)==T)

cross012$JP_NCM_A<-rowSums((cross012[c(43:48)]==0)==T)
cross012$JP_NCM_B<-rowSums((cross012[c(43:48)]==2)==T)
cross012$JP_NCM_AB<-rowSums((cross012[c(43:48)]==1)==T)

cross012$JP_NCF_A<-rowSums((cross012[c(36:41)]==0)==T)
cross012$JP_NCF_B<-rowSums((cross012[c(36:41)]==2)==T)
cross012$JP_NCF_AB<-rowSums((cross012[c(36:41)]==1)==T)

cross012$JPFP_NC_A<-rowSums((cross012[c(3:9, 11:15, 18:20, 36:48)]==0)==T)
cross012$JPFP_NC_B<-rowSums((cross012[c(3:9, 11:15, 18:20, 36:48)]==2)==T)
cross012$JPFP_NC_AB<-rowSums((cross012[c(3:9, 11:15, 18:20, 36:48)]==1)==T)

cross012$JPFP_NCM_A<-rowSums((cross012[c(3:9, 11:15, 18:20, 43:48)]==0)==T)
cross012$JPFP_NCM_B<-rowSums((cross012[c(3:9, 11:15, 18:20, 43:48)]==2)==T)
cross012$JPFP_NCM_AB<-rowSums((cross012[c(3:9, 11:15, 18:20, 43:48)]==1)==T)


cross012$JP_TXFHeterozygosity=cross012$JP_TXF_AB/(cross012$JP_TXF_AB+cross012$JP_TXF_B+cross012$JP_TXF_A)

cross012$JP_TXMHeterozygosity=cross012$JP_TXM_AB/(cross012$JP_TXM_AB+cross012$JP_TXM_B+cross012$JP_TXM_A)

cross012$JP_TXHeterozygosity=cross012$JP_TX_AB/(cross012$JP_TX_AB+cross012$JP_TX_B+cross012$JP_TX_A)



cross012$JP_NCFHeterozygosity=cross012$JP_NCF_AB/(cross012$JP_NCF_AB+cross012$JP_NCF_B+cross012$JP_NCF_A)

cross012$JP_NCMHeterozygosity=cross012$JP_NCM_AB/(cross012$JP_NCM_AB+cross012$JP_NCM_B+cross012$JP_NCM_A)

cross012$JP_NCHeterozygosity=cross012$JP_NC_AB/(cross012$JP_NC_AB+cross012$JP_NC_B+cross012$JP_NC_A)


cross012$FP_TXFHeterozygosity=cross012$FP_TXF_AB/(cross012$FP_TXF_AB+cross012$FP_TXF_B+cross012$FP_TXF_A)

cross012$FP_TXMHeterozygosity=cross012$FP_TXM_AB/(cross012$FP_TXM_AB+cross012$FP_TXM_B+cross012$FP_TXM_A)

cross012$FP_TXHeterozygosity=cross012$FP_TX_AB/(cross012$FP_TX_AB+cross012$FP_TX_B+cross012$FP_TX_A)


cross012$FP_NCFHeterozygosity=cross012$FP_NCF_AB/(cross012$FP_NCF_AB+cross012$FP_NCF_B+cross012$FP_NCF_A)

cross012$FP_NCMHeterozygosity=cross012$FP_NCM_AB/(cross012$FP_NCM_AB+cross012$FP_NCM_B+cross012$FP_NCM_A)

cross012$FP_NCHeterozygosity=cross012$FP_NC_AB/(cross012$FP_NC_AB+cross012$FP_NC_B+cross012$FP_NC_A)



cross012$JPFP_TXMHeterozygosity=cross012$JPFP_TXM_AB/(cross012$JPFP_TXM_AB+cross012$JPFP_TXM_B+cross012$JPFP_TXM_A)

cross012$JPFP_TXHeterozygosity=cross012$JPFP_TX_AB/(cross012$JPFP_TX_AB+cross012$JPFP_TX_B+cross012$JPFP_TX_A)


cross012$JPFP_NCMHeterozygosity=cross012$JPFP_NCM_AB/(cross012$JPFP_NCM_AB+cross012$JPFP_NCM_B+cross012$JPFP_NCM_A)

cross012$JPFP_NCHeterozygosity=cross012$JPFP_NC_AB/(cross012$JPFP_NC_AB+cross012$JPFP_NC_B+cross012$JPFP_NC_A)


colnames(cross012)

cross012_summary<-cross012[,c(3,5,53:132)]
head(cross012)
View(cross012_summary)
#write.csv(cross012, "11-25-2019_all_pop_data_with_AF.csv")
write.csv(cross012_summary, "12-4-2019_pop_data_summary.csv")
rm(cross012)

################ Sites expressed only in males ################

#Import .lmiss files
JoshMissingSitesNCPopFemales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesNCPopFemales.lmiss", header = T)
colnames(JoshMissingSitesNCPopFemales)<-c("CHR","POS","J_NC_F_N_DATA","J_NC_F_N_GENOTYPE_FILTERED","J_NC_F_N_MISS","J_NC_F_F_MISS")

JoshMissingSitesNCPopMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesNCPopMales.lmiss", header = T)
colnames(JoshMissingSitesNCPopMales)<-c("CHR","POS","J_NC_M_N_DATA","J_NC_M_N_GENOTYPE_FILTERED","J_NC_M_N_MISS","J_NC_M_F_MISS")

JoshMissingSitesTXPopFemales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesTXPopFemales.lmiss", header = T)
colnames(JoshMissingSitesTXPopFemales)<-c("CHR","POS","J_TX_F_N_DATA","J_TX_F_N_GENOTYPE_FILTERED","J_TX_F_N_MISS","J_TX_F_F_MISS")

JoshMissingSitesTXPopMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesTXPopMales.lmiss", header = T)
colnames(JoshMissingSitesTXPopMales)<-c("CHR","POS","J_TX_M_N_DATA","J_TX_M_N_GENOTYPE_FILTERED","J_TX_M_N_MISS","J_TX_M_F_MISS")


FelixMissingSitesXYPopMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesXYMales.lmiss", header = T)
colnames(FelixMissingSitesXYPopMales)<-c("CHR","POS","F_XY_M_N_DATA","F_XY_M_N_GENOTYPE_FILTERED","F_XY_M_N_MISS","F_XY_M_F_MISS")

FelixMissingSitesXYYPopMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/MissingSitesXYYMales.lmiss", header = T)
colnames(FelixMissingSitesXYYPopMales)<-c("CHR","POS","F_XYY_M_N_DATA","F_XYY_M_N_GENOTYPE_FILTERED","F_XYY_M_N_MISS","F_XYY_M_F_MISS")

F2MissingSitesTXMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/F2MissingSitesTXMales.lmiss", header = T)
colnames(F2MissingSitesTXMales)<-c("CHR","POS","F2_TX_M_N_DATA","F2_TX_M_N_GENOTYPE_FILTERED","F2_TX_M_N_MISS","F2_TX_M_F_MISS")
F2MissingSitesTXFemales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/F2MissingSitesTXFemales.lmiss", header = T)
colnames(F2MissingSitesTXFemales)<-c("CHR","POS","F2_TX_F_N_DATA","F2_TX_F_N_GENOTYPE_FILTERED","F2_TX_F_N_MISS","F2_TX_F_F_MISS")

F2MissingSitesNCMales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/F2MissingSitesNCMales.lmiss", header = T)
colnames(F2MissingSitesNCMales)<-c("CHR","POS","F2_NC_M_N_DATA","F2_NC_M_N_GENOTYPE_FILTERED","F2_NC_M_N_MISS","F2_NC_M_F_MISS")


F2MissingSitesNCFemales<-read.table("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/SNPs_only_expressed_M/F2MissingSitesNCFemales.lmiss", header = T)
colnames(F2MissingSitesNCFemales)<-c("CHR","POS","F2_NC_F_N_DATA","F2_NC_F_N_GENOTYPE_FILTERED","F2_NC_F_N_MISS","F2_NC_F_F_MISS")

all_missing_F2<-full_join(F2MissingSitesTXMales, F2MissingSitesTXFemales, by=c("CHR","POS"))
rm(F2MissingSitesTXMales)
rm(F2MissingSitesTXFemales)
all_missing_F2<-full_join(all_missing_F2, F2MissingSitesNCMales, by=c("CHR","POS"))
rm(F2MissingSitesNCMales)
all_missing_F2<-full_join(all_missing_F2, F2MissingSitesNCFemales, by=c("CHR","POS"))
rm(F2MissingSitesNCFemales)



#merge .lmiss files
all_missing<-full_join(JoshMissingSitesNCPopFemales, JoshMissingSitesNCPopMales, by=c("CHR","POS"))

all_missing<-full_join(all_missing, JoshMissingSitesTXPopFemales, by=c("CHR","POS"))

all_missing<-full_join(all_missing, JoshMissingSitesTXPopMales, by=c("CHR","POS"))

all_missing<-full_join(all_missing, FelixMissingSitesXYYPopMales, by=c("CHR","POS"))

all_missing<-full_join(all_missing, FelixMissingSitesXYPopMales, by=c("CHR","POS"))

rm(JoshMissingSitesNCPopFemales)
rm(JoshMissingSitesNCPopMales)
rm(JoshMissingSitesTXPopFemales)
rm(JoshMissingSitesTXPopMales)
rm(FelixMissingSitesXYYPopMales)
rm(FelixMissingSitesXYPopMales)

all_missing<-full_join(all_missing, all_missing_F2, by=c("CHR","POS"))
rm(all_missing_F2)

head(all_missing)

all_missing <- separate(all_missing, CHR, into=c("ScnbKXS", "Scaffold", "HRSCAF"), sep = "[_//;]", remove = FALSE)
colnames(all_missing)
head(all_missing)
all_missing<-all_missing[,c(3,5:45)]


all_missing_to_join<-all_missing[,c(1,2,6,10,14,18,22,26, 30, 34, 38, 42)]
head(all_missing_to_join)
write.csv(all_missing, "12-4-2019_All_missing.csv")
rm(all_missing)

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_JP_NC = case_when(
      (
        J_NC_F_F_MISS == 1 & 
          J_NC_M_F_MISS == 0) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_JP_TX = case_when(
      (
        J_TX_F_F_MISS == 1 & 
          J_TX_M_F_MISS == 0) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_FP_JP_XY = case_when(
      (
        J_TX_F_F_MISS == 1 & 
          F_XY_M_F_MISS == 0 & 
          J_TX_M_F_MISS == 0) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_FP_JP_XYY = case_when(
      (
        J_NC_F_F_MISS == 1 & 
          F_XYY_M_F_MISS == 0 & 
          J_NC_M_F_MISS == 0) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_JP_NC = case_when(
      (
        J_NC_F_F_MISS == 0 & 
          J_NC_M_F_MISS == 1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_JP_TX = case_when(
      (
        J_TX_F_F_MISS == 0 & 
          J_TX_M_F_MISS == 1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_FP_JP_XY = case_when(
      (
        J_TX_F_F_MISS == 0 & 
          F_XY_M_F_MISS == 1 & 
          J_TX_M_F_MISS == 1) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_FP_JP_XYY = case_when(
      (
        J_NC_F_F_MISS == 0 & 
          F_XYY_M_F_MISS == 1 & 
          J_NC_M_F_MISS == 1) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_JP_both = case_when(
      (
        Male_only_JP_NC == 1 & 
          Male_only_JP_TX == 1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_FP_JP_both = case_when(
      (
        Male_only_FP_JP_XY == 1 & 
          Male_only_FP_JP_XYY == 1) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_FP_JP_both = case_when(
      (
        Female_only_FP_JP_XY == 1 & 
          Female_only_FP_JP_XYY == 1) ~ 1,
      TRUE ~ 0))




colnames(all_missing_to_join)

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_F2_TX = case_when(
      (
        F2_TX_F_F_MISS > 0.9 & 
          F2_TX_M_F_MISS < 0.1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_F2_NC = case_when(
      (
        F2_NC_F_F_MISS  > 0.9  & 
          F2_NC_M_F_MISS < 0.1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_F2_both = case_when(
      (
        Male_only_F2_NC == 1 & 
          Male_only_F2_TX == 1) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_pop_and_F2_TX = case_when(
      (
        Male_only_FP_JP_XY == 1 & 
          Male_only_F2_TX == 1) ~ 1,
      TRUE ~ 0))



all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_pop_and_F2_NC = case_when(
      (
        Male_only_FP_JP_XYY == 1 & 
          Male_only_F2_NC == 1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Male_only_pop_and_F2_both = case_when(
      (
        Male_only_F2_both == 1 & 
          Male_only_FP_JP_both == 1) ~ 1,
      TRUE ~ 0))


  all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_F2_TX = case_when(
      (
        F2_TX_F_F_MISS < 0.1 & 
          F2_TX_M_F_MISS  > 0.9) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_F2_NC = case_when(
      (
        F2_NC_F_F_MISS < 0.1 & 
          F2_NC_M_F_MISS  > 0.9) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_F2_both = case_when(
      (
        Female_only_F2_NC == 1 & 
          Female_only_F2_TX == 1) ~ 1,
      TRUE ~ 0))


all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_pop_and_F2_TX = case_when(
      (
        Female_only_FP_JP_XY == 1 & 
          Female_only_F2_TX == 1) ~ 1,
      TRUE ~ 0))



all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_pop_and_F2_NC = case_when(
      (
        Female_only_FP_JP_XYY == 1 & 
          Female_only_F2_NC == 1) ~ 1,
      TRUE ~ 0))

all_missing_to_join<-   all_missing_to_join %>%
  mutate(
    Female_only_pop_and_F2_both = case_when(
      (
        Female_only_F2_both == 1 & 
          Female_only_FP_JP_both == 1) ~ 1,
      TRUE ~ 0))



head(all_missing_to_join)

colnames(all_missing_to_join)

all_missing_to_join$logical<-rowSums(all_missing_to_join[,13:35], na.rm=T)

all_missing_to_join<-subset(all_missing_to_join, all_missing_to_join$logical != 0)

rm(all_missing)
View(all_missing_to_join)
colnames(all_missing_to_join)


sum(all_missing_to_join$Male_only_JP_NC, na.rm=T)
sum(all_missing_to_join$Male_only_JP_TX, na.rm=T)
sum(all_missing_to_join$Male_only_FP_JP_XY, na.rm=T)
sum(all_missing_to_join$Male_only_FP_JP_XYY, na.rm=T)
sum(all_missing_to_join$Female_only_JP_NC, na.rm=T)
sum(all_missing_to_join$Female_only_JP_TX, na.rm=T)
sum(all_missing_to_join$Female_only_FP_JP_XY, na.rm=T)
sum(all_missing_to_join$Female_only_FP_JP_XYY, na.rm=T)
sum(all_missing_to_join$Male_only_JP_both, na.rm=T)
sum(all_missing_to_join$Male_only_FP_JP_both, na.rm=T)
sum(all_missing_to_join$Male_only_F2_TX, na.rm=T)
sum(all_missing_to_join$Male_only_F2_NC, na.rm=T)
sum(all_missing_to_join$Male_only_F2_both, na.rm=T)
sum(all_missing_to_join$Male_only_pop_and_F2_TX, na.rm=T)
sum(all_missing_to_join$Male_only_pop_and_F2_NC, na.rm=T)
sum(all_missing_to_join$Male_only_pop_and_F2_both, na.rm=T)
sum(all_missing_to_join$Female_only_F2_TX, na.rm=T)
sum(all_missing_to_join$Female_only_F2_NC, na.rm=T)
sum(all_missing_to_join$Female_only_F2_both, na.rm=T)
sum(all_missing_to_join$Female_only_pop_and_F2_TX, na.rm=T)
sum(all_missing_to_join$Female_only_pop_and_F2_NC, na.rm=T)
sum(all_missing_to_join$Female_only_FP_JP_both, na.rm=T)
sum(all_missing_to_join$Female_only_pop_and_F2_both, na.rm=T)



write.csv(all_missing_to_join, "12-3-2019_missing_data.csv")
rm(all_missing_to_join)

###########################








#####*************** Join datasets ***************#####

crossdata<-read.csv("12-3-2019_pedigree_summary.csv")
popdata<-read.csv("12-4-2019_pop_data_summary.csv")
F2data<-read.csv("12-3-2019_F2_summary.csv")
male_only_sites<-read.csv("12-3-2019_missing_data.csv")

rm(crossdata)
rm(F2data)
rm(popdata)
head(crossdata)
head(popdata)
head(F2data)


#crossdata<-crossdata[-1,]
#popdata<-popdata[-1,]


joined_data<-full_join(crossdata, F2data, by=c("Scaffold", "Position"))

joined_data<-full_join(joined_data, popdata, by=c("Scaffold", "Position"="POS"))

joined_data<-full_join(joined_data, male_only_sites, by=c("Scaffold", "Position"="POS"))

head(joined_data)
write.csv(joined_data, "12-4-19-All_SNPs_joined.csv")





#####********** Convert all scaffold positions to chromosome positions using Python script SNP_list_Interactive_Chromonomer_LG_position_switcher.py **********#####
############# TX ##################
Joined_dataset
Joined_dataset<-Joined_dataset%>% mutate(Position_adjusted=Position)
head(Joined_dataset)
#Joined_dataset<-left_join(Joined_dataset, TX_LGs, "Scaffold", Fill=T)
str(Joined_dataset)
str(TX_LGs)
##View(Joined_dataset)

############ Join finalized map to reordered markers ######################
#TX_stripped_order<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/Sex_linked_SNPs/TX_stripped_order.csv")
#TX_stripped_order<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/Sex_linked_SNPs/TX_stripped_order.csv")
#converted_test<-read.csv("testSNPsconverted.csv", header = T, stringsAsFactors = F)
#TX_stripped_order<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/Sex_linked_SNPs/TX_stripped_order.csv")
TX_final_order<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/ASMAP_final/Map_manual_edits/Final_edited_markers_positions_colocated_filled.csv")
TX_final_order<-read.csv("Final_edited_markers_positions_colocated_filled.csv")


Joined_dataset_LGs<-left_join(converted_test,TX_final_order, by=c("Scaffold"="scaffold","Position"="position"), fill=T)


Joined_dataset_LGs<-left_join(Joined_dataset,TX_stripped_order, by=c("Scaffold"="CHROM"), fill=T)


#colnames(Joined_dataset)<-c("X","CHROM", "Position","LG","CM" , "Position_adjusted")

#View(Joined_dataset_LGs)
Joined_dataset_LGs<-left_join(Joined_dataset,TX_stripped_order, by=c("Scaffold"="CHROM"), fill=T)
"D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/mapped_and_colocated_markers_TX.csv")
linkage_map_markers<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/mapped_and_colocated_markers_TX.csv")
linkage_map_markers<-read.csv("D:/Dropbox/Dropbox/Professional/University_of_Toronto/Genomics/HiCSNPs/TranscriptomeLinkageMap/ASMAP/TX_transcriptome/mapped_and_colocated_markers_TX.csv")

head(linkage_map_markers)


Joined_dataset_LGs<-left_join(Joined_dataset_LGs,linkage_map_markers, by = c("Scaffold","Position"))

######################## Summarize data ################################

#Joined_dataset<-read.csv("line2removed12-11-19-All_SNPs_joined.csv")

#Converted_dataset<-read.csv("12-17-19-All_SNPs_joined_converted.csv", header = T, stringsAsFactors = F)

#converted_test<-read.csv("testSNPsconverted.csv", header = T, stringsAsFactors = F)


#converted_test %>% group_by(LG) %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 

#summary_all_columns<-Converted_dataset %>% group_by(LG) %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(summary_all_columns, "Summary_all_SNPs_by_LG_12-17-2019.csv")


#colSums(Joined_dataset, na.rm=T)


#zeros<-Converted_dataset %>% group_by(LG) %>% filter(LG_position ==0) 
#View(zeros)
Windowed_summary_lengths10K <- Converted_dataset %>% group_by(LG)  %>% mutate(position_window=LG_position%/%10000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(Windowed_summary_lengths10K, "10KWindows.csv")

#windoweddata<-read.csv("10KWindows.csv")
#windoweddata<-data.frame(windoweddata)

#test<-subset(Windowed_summary_lengths10K, Windowed_summary_lengths10K$LG=="L.10")
#plot(test$LG_position_mean, test$n_sum)

#plot(test$LG_position_mean, test$JXTX_YLinked_sum)
#plot(test$LG_position_mean, test$JXNC_Autosomal_sum)
#plot(test$LG_position_mean, test$F2_XY_pattern.y_sum)
#plot(test$LG_position_mean, test$JP_TXMHeterozygosity_mean)
#plot(test$LG_position_mean, (test$JP_TXMHeterozygosity_mean-test$JP_TXFHeterozygosity_mean))

#hist(test$Position_mean)

Windowed_summary_lengths50K <- Converted_dataset %>% group_by(LG) %>% mutate(position_window=LG_position%/%50000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(Windowed_summary_lengths50K, "50KWindows.csv")

Windowed_summary_lengths100K <- Converted_dataset %>% group_by(LG) %>% mutate(position_window=LG_position%/%100000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(Windowed_summary_lengths100K, "100KWindows.csv")

Windowed_summary_lengths200K <- Converted_dataset %>% group_by(LG) %>% mutate(position_window=LG_position%/%200000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(Windowed_summary_lengths200K, "200KWindows.csv")


Windowed_summary_lengths500K <- Converted_dataset %>% group_by(LG) %>% mutate(position_window=LG_position%/%500000) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs(sum(., na.rm = T), mean(., na.rm = T))) 
write.csv(Windowed_summary_lengths500K, "500KWindows.csv")



#rm(cross012)
joshcrossdata<-cross012
nrow(joshcrossdata)
nrow(subset(joshcrossdata, is.na(joshcrossdata$sampleTXDAD_)==F))
#View(joshcrossdata)
sum(is.na(joshcrossdata))
length(unique(joshcrossdata$CHROM))
sum(Joined_dataset_LGs$JXTX_YLinked, na.rm=T)
sum(Joined_dataset_LGs$JXTX_XLinked, na.rm=T)
sum(Joined_dataset_LGs$JXTX_Hemizygous, na.rm=T)
sum(Joined_dataset_LGs$JXTX_Autosomal, na.rm=T)

sum(Joined_dataset_LGs$JXNC_YLinked, na.rm=T)
sum(Joined_dataset_LGs$JXNC_XLinked, na.rm=T)
sum(Joined_dataset_LGs$JXNC_Hemizygous, na.rm=T)
sum(Joined_dataset_LGs$JXNC_Autosomal, na.rm=T)

sum(Joined_dataset_LGs$JX_shared_YLinked, na.rm=T)
sum(Joined_dataset_LGs$JX_shared_XLinked, na.rm=T)
sum(Joined_dataset_LGs$JX_shared_Hemizygous, na.rm=T)
sum(Joined_dataset_LGs$JX_shared_Autosomal, na.rm=T)


colnames(Joined_dataset_LGs)
Joined_dataset_LGs_in_map<-subset(Joined_dataset_LGs, is.na(Joined_dataset_LGs$Position_adjusted_LG)==F)
#View(unique(Joined_dataset_LGs_in_map$CHROM))
sum(Joined_dataset_LGs_in_map$JXTX_YLinked, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXTX_XLinked, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXTX_Hemizygous, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXTX_Autosomal, na.rm=T)

sum(Joined_dataset_LGs_in_map$JXNC_YLinked, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXNC_XLinked, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXNC_Hemizygous, na.rm=T)
sum(Joined_dataset_LGs_in_map$JXNC_Autosomal, na.rm=T)


sum(Joined_dataset_LGs$JP_male_only, na.rm=T)
sum(Joined_dataset_LGs_in_map$JP_male_only, na.rm=T)

########### Windowed recombination rates #######################
library(dplyr)
options(scipen=999)

setwd("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Windowed_analyses/")
final_map_converted<-read.table("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Chromonomer_final/CHRR_linkage_map.tsv", sep="\t", stringsAsFactors = F)
colnames(final_map_converted)<-c("LG", "Scaff_full", "CM_position","broken_contig","orientation","broken_contig_position","Position_on_LG","broken_contig_ordinal") ## TX version


final_map_converted<-read.csv("E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/Rumex_genome_paper/Finalized_analyses/Marey_Maps/AGP_converter/TX_AGP_new_final_edited_NC_old_map1217.csv", stringsAsFactors = F, header = T)
colnames(final_map_converted)<-c("scaffold","scaffold_position","LG", "Position_on_LG", "Old_LG", "CM_position") ## NC version

head(final_map_converted)
getwd()

final_map_converted<-select(final_map_converted, c(LG,Position_on_LG,CM_position))
#plot(subset(final_map_converted, final_map_converted$LG=="L.8")$Position_on_LG, subset(final_map_converted, final_map_converted$LG=="L.8")$CM, main = "Marey Map, A2", ylab = "CM_position", xlab="Start position along LG")

window_size<-50000
window_size<-100000
window_size<-500000
window_size<-1000000
window_size<-5000000

assign(paste("finalmap_",as.character(window_size),"_bp_recombination_rates", sep=""), final_map_converted %>% group_by(LG) %>% mutate(position_window=Position_on_LG%/%as.numeric(window_size)) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% do (model = lm(CM_position~Position_on_LG, data = .)) %>% mutate(coefficient = coef(model)[2]))
ls()

finalmap_windows <- final_map_converted %>% group_by(LG) %>% mutate(position_window=Position_on_LG%/%window_size) %>% group_by(LG,position_window) %>% add_tally() %>% select_if(., is.numeric) %>% summarize_all(funs( mean(., na.rm = T)))

assign(paste("full_map_",as.character(window_size),"_bp_recombination_rates", sep=""), left_join(finalmap_windows, (eval(parse(text = paste("finalmap_",as.character(window_size),"_bp_recombination_rates", sep="") ))), by=c("LG","position_window")))

#View(current_map)
current_map<-(eval(parse(text = paste("full_map_",as.character(window_size),"_bp_recombination_rates", sep=""))))  
       
#str(current_map)
colnames(current_map)
current_map<-current_map[,c(1:5,7)]
#current_map<-current_map[,c(1:3,5)]


#write.csv(current_map, paste("TX_full_map_",as.character(window_size),"_bp_recombination_rates_",Sys.Date(), ".csv", sep=""))
#write.csv(current_map, paste("NC_full_map_",as.character(window_size),"_bp_recombination_rates_",Sys.Date(), ".csv", sep=""))

plot(subset(current_map, current_map$LG=="L.10")$Position_on_LG, subset(current_map, current_map$LG=="L.10")$n)
plot(subset(current_map, current_map$LG=="L.10")$Position_on_LG, subset(current_map, current_map$LG=="L.10")$coefficient)

plot(subset(current_map, current_map$LG=="L.10")$Position_on_LG, subset(current_map, current_map$LG=="L.10")$n)
plot(subset(current_map, current_map$LG=="L.10")$Position_on_LG, subset(current_map, current_map$LG=="L.10")$coefficient, main = "Recombination by 5MB window, X", ylab = "Coefficient of linear model", xlab="Start position along LG")


plot(subset(current_map, current_map$LG=="L.8")$Position_on_LG, subset(current_map, current_map$LG=="L.8")$n)
plot(subset(current_map, current_map$LG=="L.8")$Position_on_LG, subset(current_map, current_map$LG=="L.8")$coefficient, main = "Recombination by 5MB window, A2", ylab = "Coefficient of linear model", xlab="Start position along LG")


SNP_count<-10
SNP_count<-50
SNP_count<-100

assign(paste("finalmap_",as.character(SNP_count),"_snp_recombination_rates", sep=""),final_map_converted %>% group_by(LG) %>%   select_if(., is.numeric) %>% arrange(., Position_on_LG) %>% mutate(obs_count = row_number()) %>% mutate(snp_count_window=obs_count%/%SNP_count) %>% group_by(LG, snp_count_window ) %>% do (model = lm(CM_position~Position_on_LG, data = .)) %>% mutate(coefficient = coef(model)[2]))

finalmap_snp_windows<- final_map_converted %>% group_by(LG) %>%   select_if(., is.numeric) %>% arrange(., Position_on_LG) %>% mutate(obs_count = row_number()) %>% mutate(snp_count_window=obs_count%/%SNP_count) %>% group_by(LG, snp_count_window )  %>% summarise(., mean(Position_on_LG))

#View(finalmap_snp_windows)
assign(paste("full_map_",as.character(SNP_count),"_snp_windows", sep=""),left_join(finalmap_snp_windows, (eval(parse(text = paste( "finalmap_",as.character(SNP_count),"_snp_recombination_rates", sep="")))) , by=c("LG","snp_count_window")))


current_map<-(eval(parse(text = paste("full_map_",as.character(SNP_count),"_snp_windows", sep=""))))  
       
colnames(current_map)
#current_map<-current_map[,c(1:5,7)]
current_map<-current_map[,c(1:3,5)]
#View(current_map)
#write.csv(current_map, paste("TX_full_map_",as.character(SNP_count),"_snp_windows_",Sys.Date(), ".csv", sep=""))
#write.csv(current_map, paste("NC_full_map_",as.character(SNP_count),"_snp_windows_",Sys.Date(), ".csv", sep=""))


plot(subset(current_map, current_map$LG=="L.8")$Position_on_LG, subset(current_map, current_map$LG=="L.8")$coefficient, main = "Recombination by 10-snp window, A2", ylab = "Coefficient of linear model", xlab="Start position along LG")

plot(subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.10")$Position_on_LG, subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.10")$coefficient, main = "Recombination by 50-snp window, X", ylab = "Coefficient of linear model", xlab="Start position along LG")

plot(subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.3")$Position_on_LG, subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.3")$coefficient)

plot(subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.5")$Position_on_LG, subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.5")$coefficient)

plot(subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.7")$Position_on_LG, subset(finalmap_50_SNP_windows, finalmap_50_SNP_windows$LG=="L.7")$coefficient)


