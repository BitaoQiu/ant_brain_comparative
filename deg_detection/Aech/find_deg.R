library(base)
library(tximport)
library('DESeq2')
#aech_files = c(paste('Aech',c(1,2,'3x',4:9,'10x',11,'12x','13x',14,15),sep ='_'))
aech_files = c(paste('Aech',c(4:9,'10x',11,'12x','13x',14,15),sep ='_'))
files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Aech_gemoma_t2g.txt',header = F, sep = '\t',stringsAsFactors = F)
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
cor(txi.salmon$abundance,method = 'p')[seq(3,12,3),seq(3,12,3)]
batch = factor(c(rep(c('gyne','major','minor'),4)))#,     c('minor','gyne','minor','gyne')))
colony = factor(c(rep(c(1:4),each = 3)))#,c(1,4,5,5)))
sampleTable <- data.frame(condition = batch, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition  + colony)# + type )
dds <- DESeq(dds)
res_gb <- results(dds, contrast=c("condition","gyne","major"),alpha = 0.05)
res_gb <- lfcShrink(dds, contrast=c("condition","gyne","major"), res=res_gb)
res_gs <- results(dds, contrast=c("condition","gyne","minor"),alpha = 0.05)
res_gs = lfcShrink(dds, contrast=c("condition","gyne","minor"), res=res_gs)
write.table(res_gb, "Aech_gyne_major.csv",quote = F, row.names = T, col.names = T, sep = '\t')
write.table(res_gs, "Aech_gyne_minor.csv",quote = F, row.names = T, col.names = T, sep = '\t')
summary(res_gb)
res_bs <- results(dds, contrast=c("condition","major","minor"),alpha = 0.05)
cut_off = 0.05
gyne_bias = intersect(rownames(res_gb[which(res_gb$padj < cut_off  & res_gb$log2FoldChange > 0),]),rownames(res_gb[which(res_gs$padj < cut_off & res_gs$log2FoldChange > 0),]))
major_bias = intersect(rownames(res_gb[which(res_gb$padj < cut_off & res_gb$log2FoldChange < 0),]),rownames(res_bs[which(res_bs$padj < cut_off & res_bs$log2FoldChange > 0),]))
minor_bias = intersect(rownames(res_gs[which(res_gs$padj < cut_off & res_gs$log2FoldChange < 0),]),rownames(res_bs[which(res_bs$padj < cut_off & res_bs$log2FoldChange < 0),]))
worker_bias = union(rownames(res_gb[which(res_gb$padj < cut_off & res_gb$log2FoldChange < 0),]),rownames(res_gs[which(res_gs$padj < cut_off & res_gs$log2FoldChange < 0),]))
worker_bias = worker_bias[!worker_bias %in% c(major_bias,minor_bias)]
#####
gyne_bias = rownames(res_gs[which(res_gs$padj < 0.1 & res_gs$log2FoldChange > 0),])
worker_bias = rownames(res_gs[which(res_gs$padj < 0.1 & res_gs$log2FoldChange < 0),])
#####
gyne_bias = rownames(res_gb[which(res_gb$padj < 0.1 & res_gb$log2FoldChange > 0),])
worker_bias = rownames(res_gb[which(res_gb$padj < 0.1 & res_gb$log2FoldChange < 0),])
####
major_bias_2 = rownames(res_bs[which(res_bs$padj < 0.1 & res_bs$log2FoldChange > 0),])
minor_bias_2 = rownames(res_bs[which(res_bs$padj < 0.1 & res_bs$log2FoldChange < 0),])

