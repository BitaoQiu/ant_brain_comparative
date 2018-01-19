library(tximport)
library('DESeq2')
aech_files = c(paste('Aech',c(1:15),sep ='_'), paste('Aech',c('3x','10x','12x','13x'),sep = '_'))
aech_files = c(paste('Aech',c(1:2,'3x',4:9,'10x',11,'12x','13x',14,15),sep ='_'))

files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Aech_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
plot(hclust(dist(t((txi.salmon$abundance +1 )))))
abundance = txi.salmon$abundance
gyne_major = c()
gyne_minor = c()
major_minor = c()
for(i in seq(1,5,1)){
  index = (i-1)*3 + 1
  gyne_major[i] = cor(abundance[,index],abundance[,(index+1)],method = 'k')
  gyne_minor[i] = cor(abundance[,index],abundance[,(index+2)],method = 'k')
  major_minor[i] = cor(abundance[,(index+1)],abundance[,(index+2)],method = 'k')}
gyne_exp = matrix(abundance[,seq(1,15,3)],ncol = 1)
major_exp = matrix(abundance[,seq(2,15,3)],ncol = 1)
minor_exp = matrix(abundance[,seq(3,15,3)], ncol = 1)
cor(gyne_exp,major_exp,method = 'k')
cor(gyne_exp,minor_exp,method = 'k')
cor(major_exp,minor_exp,method = 'k')


# Between gyne and workers:
batch = factor(c(rep(c('gyne','major','minor'),5)))#,     c('minor','gyne','minor','gyne')))
colony = factor(c(rep(c(1:5),each = 3)))#,c(1,4,5,5)))
sampleTable <- data.frame(condition = batch, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition  + colony)# + type )
dds <- dds[ rowSums(counts(dds)) > 15, ]
dds <- DESeq(dds)
######
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:4000]

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
df <- as.data.frame(colData(dds)[,c("condition",  'colony')]) 
select <- sample(dim(dds)[1],4000)
pdf( "Aech_cluster.pdf", width = 5, height = 4 )
pheatmap(assay(vsd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T,annotation_col = sampleTable)
plotPCA(vsd, intgroup=c('condition'), ntop = 200)

res <- results(dds, contrast=c("condition","major","minor"))
plotMA(res, ylim = c(-4,4))
summary(res)
table(abs(res$log2FoldChange) > 1)
length(which(res$padj< 0.1))
summary(res$pad)
#####

res[order((res$log2FoldChange),decreasing = T)[c(1:10)],]
test_gene = 'gene_4357'
plotCounts(dds, gene=test_gene, intgroup="condition",col = batch, pch = 19)
res[test_gene,]
tx2gene[which(tx2gene$V2== test_gene),]
#########