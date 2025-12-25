library(DESeq2)

df <- read.table("../GSE131592_counts.txt")
df <- df[,-1:-5]

#make metadata
meta_df <- data.frame(sample=colnames(df),
                      condition=factor(),
                      row.names=colnames(df))
meta_df$sex = "M"
meta_df$sex[df["Xist",]>10] <- "F"
meta_df$sex = factor(meta_df$sex,levels=c("M","F"))

#set up object
df <- sapply(df,as.integer)
dds <- DESeqDataSetFromMatrix(countData=df,colData=meta_df,design=~condition+sex)

#remove useless features
dds <- dds[rowSums(counts(dds)>2)>2,]

#DESeq
dds <- DESeq(dds)

res <- as.data.frame(results(dds, contrast = c("condition","KC","KPC"),cooksCutoff=FALSE, independentFiltering=FALSE))
res1a <- as.data.frame(lfcShrink(dds,coef=length(resultsNames(dds)),
                                 res=results(dds, contrast = c("condition","KC","KPC"),cooksCutoff=FALSE,
                                             independentFiltering=FALSE),type = "ashr"))

write.table(res1,"DE.genes.txt")
write.table(res1a,"shDE.genes.txt")
