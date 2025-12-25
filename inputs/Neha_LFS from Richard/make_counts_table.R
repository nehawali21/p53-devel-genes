proj <- c("GSE131592","GSE162698","GSE163088","GSE226389","GSE226406","GSE36952")

for(t in 1:length(proj)){
  
  fm <- paste0("/mnt/lustre/processed/Lu_lab/RichardOwen/RNASeq/Results/",proj[t],"/1.4.1/featureCounts/")
  fm <- list.files(fm,full.names = TRUE)
  fm <- fm[!grepl(".summary",fm)]
  
  x <- NULL
  df <- read.table(fm[1])
  
  for(i in 2:length(fm)){
    x <- read.table(fm[i])
    df <- cbind(df,x[,7])
    
  }
  
  colnames(df) <- df[1,]
  df <- df[-1,]
  rownames(df) <- df[,1]
  df <- df[,-1]
  
  colnames(df) <- substr(colnames(df),1,10)
  
  write.table(df,paste0(proj[t],"_counts.txt"),quote = FALSE)
  
}

#make samples table
df <- NULL

for(i in 1:length(proj)){
  x <- read.table(paste0("/mnt/lustre/processed/Lu_lab/RichardOwen/RNASeq/Results/",proj[i],"/1.4.1/sample_names.tsv"))
  x$project <- proj[i]
  df <- rbind(df,x)
}

colnames(df) <- df[1,]
df <- df[-1,
]
write.table(df,"sample_names.txt",quote=FALSE)