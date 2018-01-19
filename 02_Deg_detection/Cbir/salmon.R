library(tximport)
library('DESeq2')
lhum_files = paste('Cbir',c(1:8),sep ='_')

files <- file.path('quants', lhum_files, "quant.sf")
names(files) <- lhum_files
tx2gene <- read.csv('Cbir_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
# Between gyne and workers:
batch = factor(rep(c('gyne','worker'),4))
colony = factor(rep(c(1:4),each = 2))
sampleTable <- data.frame(condition = batch, colony = colony)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition + colony)
#dds <- dds[ rowSums(counts(dds)) > 8, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','gyne','worker'),alpha = 0.05)
res <- lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
summary(res)
write.table(res, "Cbir_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')
summary(res)
######
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]
vsd <- varianceStabilizingTransformation(dds, blind=F)
df <- as.data.frame(colData(dds)[,c("condition",'colony')]) 
select <- sample(dim(dds)[1],4000)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col=sampleTable)
plotPCA(vsd, intgroup=c("condition"),ntop = 200)

plotMA(res)
table(abs(res$log2FoldChange) > 1)
summary(res$pad<0.1)
summary(res$pad)
summary(res)
#####

res[order(abs(res$log2FoldChange),decreasing = T)[c(1:10)],]
test_gene = 'gene_14666'
plotCounts(dds, gene=test_gene, intgroup="condition",col = batch, pch = 19)
res[test_gene,]
#########
gyne_bias = rownames(res[which(res$padj < 0.1 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.1 & res$log2FoldChange < 0),])

