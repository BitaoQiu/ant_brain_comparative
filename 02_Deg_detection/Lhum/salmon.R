library(tximport)
library('DESeq2')
lhum_files = paste('Lhum',c(1:10),sep ='_')

files <- file.path('quants', lhum_files, "quant.sf")
names(files) <- lhum_files
tx2gene <- read.csv('Lhum_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
plot(hclust(dist(t(log(txi.salmon$counts +1)))))
abundance = txi.salmon$abundance

# Between gyne and workers:
batch = factor(rep(c('gyne','worker'),5))
colony = factor(rep(c(1:5),each = 2))
sampleTable <- data.frame(condition = batch, colony = colony)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition + colony)
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','gyne','worker'))
######
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("condition",'colony')]) 
select <- sample(dim(dds)[1],4000)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col=sampleTable)
plotPCA(vsd, intgroup=c("condition"),ntop = 200)

plotMA(res)
summary(res)
table(abs(res$log2FoldChange) > 1)
summary(res$pad<0.1)
summary(res$pad)
#####

res[order((res$log2FoldChange),decreasing = F)[c(1:10)],]
test_gene = 'gene_5081'
plotCounts(dds, gene=test_gene, intgroup="condition",col = batch, pch = 19)
res[test_gene,]
tx2gene[tx2gene$V2 == test_gene,]
#########