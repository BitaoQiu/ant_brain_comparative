library(base)
library(tximport)
library('DESeq2')
aech_files = c(paste('Lhum',c(1:10),sep ='_'))
files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Lhum_gemoma_t2g.txt',header = F, sep = '\t',stringsAsFactors = F)
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
batch = factor(c(rep(c('gyne','worker'),4)))#,     c('minor','gyne','minor','gyne')))
colony = factor(c(rep(c(2:5),each = 2)))#,c(1,4,5,5)))
batch = factor(c(rep(c('gyne','worker'),5)))#,     c('minor','gyne','minor','gyne')))
colony = factor(c(rep(c(1:5),each = 2)))#,c(1,4,5,5)))

sampleTable <- data.frame(condition = batch, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition  + colony)# + type )
dds = DESeq(dds)
res <- results(dds,contrast = c('condition','gyne','worker'),alpha = 0.05)
summary(res)
res <- lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
write.table(res, "Lhum_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')

summary(res)
#res = lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
gyne_bias = rownames(res[which(res$padj <  0.05 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.05 & res$log2FoldChange < 0),])

