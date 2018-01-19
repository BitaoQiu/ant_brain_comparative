library(base)
library(tximport)
library('DESeq2')
aech_files = c(paste('Pcan',c('134','65b','73b','75b','76b','127','137','135','66b','72b'),sep ='_'))
files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Pcan_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
caste = factor(c(rep('gyne',4),rep('worker',6)))#,     c('minor','gyne','minor','gyne')))
mated = factor(c(rep('a',4),'b',rep('a',2),rep('b',3)))
sampleTable <- data.frame(caste = caste, mated = mated)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~caste + mated)# + type )
dds <- DESeq(dds)
res <- results(dds,contrast = c('caste','gyne','worker'))

res <- lfcShrink(dds, contrast=c("caste","gyne","worker"), res=res)
plotMA(res)
write.table(res, "Pcan_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')

gyne_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange < 0),])
summary(res)
ant_res = read.table('ant_res.txt',header = T,sep = '\t')
gene_ortholog_table = read.table('../../ortholog/ortholog_ant_bee/gene_table.proteinortho', col.names = c('Aech','Sinv','Mpha','Lhum','Cbir','Dqua','Amel'))
#gene_ortholog_table$Amel = paste('Amel',gene_ortholog_table$Amel,sep = '_')
res.data = data.frame(res)
res.data$Aech_ID = as.character(gene_ortholog_table$Aech)[match(row.names(res.data),gene_ortholog_table$Amel)]

#res = data.frame(res)
table(rownames(ant_res[which(ant_res$padj < 0.005),]) %in% rownames(res[which(res[,5] < 0.005),]))
res[rownames(ant_res[which(ant_res$log2FoldChange > 0 & rownames(ant_res) %in% rownames(res[which(res[,5] < 0.01 & res[,2] > 0),])),]),]
ant_res[which(ant_res$log2FoldChange < 0 & rownames(ant_res) %in% rownames(res[which(res[,5] < 0.01 & res[,2] < 0),])),]
ant_res[which(ant_res$log2FoldChange > 0 & rownames(ant_res) %in% rownames(res[which(res[,5] < 0.0001 & res[,2] > 0),])),]
dim(ant_res[which(ant_res$log2FoldChange > 0 & rownames(ant_res) %in% rownames(res[which(res[,5] < 0.0001 & res[,2] > 0),])),])
dim(ant_res[which(ant_res$log2FoldChange < 0 & rownames(ant_res) %in% rownames(res[which(res[,5] < 0.0001 & res[,2] < 0),])),])

ant_res$log2_bee = res.data$log2FoldChange[match(rownames(ant_res),res.data$Aech_ID)]
ant_res$padj_bee = res.data$padj[match(rownames(ant_res),res.data$Aech_ID)]
tmp = ant_res[,c(7,2,8,6,9)]
names(tmp) = c("Gene_name",'Log2_ant','Log2_bee','padj_ant','padj_bee')
tmp = tmp[order(tmp$Log2_ant*tmp$Log2_bee,decreasing = T),]
tmp
write.table(tmp,'ant_bee_deg.tsv',quote = F,sep = '\t',row.names = F, col.names = T)
