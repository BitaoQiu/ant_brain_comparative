library(tximport)
library('DESeq2')
lhum_files = paste('Dqua',c(1:5,7:13),sep ='_')

files <- file.path('quants', lhum_files, "quant.sf")
names(files) <- lhum_files
tx2gene <- read.csv('Dqua_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
# Between gyne and workers:
batch = factor(c(rep('gyne',6),rep('worker',6)))
colony = factor(rep(c(24:29),2))
sampleTable <- data.frame(condition = batch, colony = colony)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition + colony)
#dds <- dds[ rowSums(counts(dds)) > 13, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','gyne','worker'),alpha = 0.05)
summary(res)
res <- lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
write.table(res, "Dqua_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')

summary(res)
######
ant_bee = read.table('../Amel/ant_bee_deg.tsv',header = T,sep = '\t')
annotate = read.table('../../ortholog/gene_table.poff')
res.data = data.frame(res)
rownames(res.data) = paste('Dqua',rownames(res.data) ,sep = '_')
res.data$Aech_ID = annotate$V1[match(rownames(res.data),annotate$V6)]
ant_bee$Log2_Dqua = res.data$log2FoldChange[match(row.names(ant_bee),res.data$Aech_ID)]
ant_bee$padj_Dqua = res.data$padj[match(row.names(ant_bee),res.data$Aech_ID)]
ant_bee[order(ant_bee$Log2_Dqua*ant_bee$Log2_ant,decreasing = T),c(1,2,3,6,7)]
######
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("condition")])#,'colony')]) 
select <- sample(dim(dds)[1],4000)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col=sampleTable)
plotPCA(vsd, intgroup=c("condition"),ntop = 200)

plotMA(res, ylim = c(-4,4))
table(abs(res$log2FoldChange) > 1)
summary(res$pad<0.1)
summary(res$pad)
summary(res)
#####

res[order(abs(res$log2FoldChange),decreasing = T)[c(1:10)],]
test_gene = 'gene_11190'
plotCounts(dds, gene=test_gene, intgroup="condition",col = batch, pch = 19)
res[test_gene,]
tx2gene[tx2gene$V2 == test_gene,]
#########
gyne_bias = rownames(res[which(res$padj < 0.1 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.1 & res$log2FoldChange < 0),])
