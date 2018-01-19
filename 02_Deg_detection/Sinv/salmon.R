library(tximport)
library('DESeq2')
lhum_files = paste('Sinv',c(1:16),sep ='_')

files <- file.path('quants', lhum_files, "quant.sf")
names(files) <- lhum_files
tx2gene <- read.csv('Sinv_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
plot(hclust(dist(t(log(txi.salmon$counts +1)))))
abundance = txi.salmon$abundance
gyne_worker = c()
for(i in seq(1,4,1)){
  index = (i-1)*4 + 1
  gyne_worker[i] = cor(abundance[,(index+2)],abundance[,(index+3)],method = 'k')}

gyne_exp = matrix(abundance[,seq(3,16,4)], ncol = 1)
worker_exp = matrix(abundance[,seq(4,16,4)], ncol = 1)
cor(gyne_exp,worker_exp,method = 'k')
# Between gyne and workers:
batch = factor(rep(c('PVQ','PW','AVQ','AW'),4))
colony = factor(rep(c(1:4),each = 4))
sampleTable <- data.frame(condition = batch, colony = colony)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition + colony)
dds <- dds[ rowSums(counts(dds)) > 16, ]
dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','AVQ','AW'))
######
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("condition",'colony')]) 
select <- sample(dim(dds)[1],4000)
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col=sampleTable)
plotPCA(vsd, intgroup=c("condition"),ntop = 15000)
res <- results(dds,contrast = c('condition','AVQ','AW'))
plotMA(res)
summary(res)
table(abs(res$log2FoldChange) > 1)
summary(res$pad<0.1)
summary(res$pad)
#####

res[order(abs(res$log2FoldChange),decreasing = T)[c(1:10)],]
test_gene = 'gene_13551'
plotCounts(dds, gene=test_gene, intgroup="condition",col = batch, pch = 19)
res[test_gene,]
tx2gene[tx2gene$V2 == test_gene,]
#########