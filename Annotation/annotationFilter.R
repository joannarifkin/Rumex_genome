
base <- fread('base.txt',header = FALSE)
base <- data.frame(base$V1[-length(base$V1)])
names(base) <- c("base") 

reads <- base
dna <- c("TXOAK20m","TXROS2f","SCMAR2f","SCMAR1m")

loops <- 1
for (ind in dna){
  cat(ind,"\n")
  ind <- toString(ind)
  fileName <- paste(ind,".gen.reads.txt",sep="")
  file <- fread(fileName,header = FALSE)
  file$cov <-  (file$V3/(file$V2/1000)) / sum(file$V3/1000000) 
  
  file <- file[-length(file$V1),]
  file <- file[,c(1,5)]
  names(file) <- c("locus",ind)
  reads <- data.frame(sqldf('select reads.*, file.* from reads left join file on reads.base = file.locus'))
  columns <- c(1:loops,(loops+2))
  if (loops == 1){
    columns <- c(1,3)
  }
  loops <- loops + 1
  reads <- reads[,columns]
  #print(length(reads$V1))
}

reads <- reads[apply(reads[,-1], 1, function(x) all(x<(-log(0.05)))),]
genreads <- reads
#reads <- reads[reads$TXOAK20m < -log(0.05),]
#ggplot(reads,aes(x=TXOAK20m)) + geom_histogram() 
#genreads$TXOAK20m[genreads$base == "IDg66968"]

reads <- data.frame(reads[,1])
names(reads) <- c("base") 

rna <- c("pollen_TX1B","pollen_TX2B","flower_17_TXMale","flower_24_TXFem","leaf_TXROS24M","leaf_TX1F")

loops <- 1
for (ind in rna){
  cat(ind,"\n")
  ind <- toString(ind)
  fileName <- paste(ind,".rna.reads.txt",sep="")
  file <- fread(fileName,header = FALSE)
  file$cov <-  (file$V3/(file$V2/1000)) / sum(file$V3/1000000) 
  
  file <- file[-length(file$V1),]
  file <- file[,c(1,5)]
  names(file) <- c("locus",ind)
  reads <- data.frame(sqldf('select reads.*, file.* from reads left join file on reads.base = file.locus'))
  columns <- c(1:loops,(loops+2))
  if (loops == 1){
    columns <- c(1,3)
  }
  loops <- loops + 1
  reads <- reads[,columns]
  #print(length(reads$V1))
}

#allExp <- reads[apply(reads[,-1], 1, function(x) !all(x==0)),]
allExp <- reads[apply(reads[,-1], 1, function(x) any(x>=0.03)),]

ggplot(allExp,aes(x=pollen_TX1B)) + geom_histogram() 

write.table(allExp$base, file = "goodgenes.list", append = FALSE, quote = FALSE, sep = "\t", 
            eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
            col.names = FALSE, qmethod = c("escape", "double"))

ggplot()
