library(base)
library(tximport)
library('DESeq2')
aech_files = c(paste('Amel',c(1:6),sep ='_'))
files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Amel.t2g.txt',header = F, sep = '\t',skip = 1)
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
batch = factor(c(rep(c('worker','gyne'),each =3 )))#,     c('minor','gyne','minor','gyne')))

sampleTable <- data.frame(condition = batch)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition)# + type )
dds <- DESeq(dds)

res <- results(dds, c('condition','gyne','worker'))
#res = lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
gyne_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange < 0),])
summary(res)

