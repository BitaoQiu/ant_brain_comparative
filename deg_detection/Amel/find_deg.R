library(base)
library(tximport)
library('DESeq2')
aech_files = c(paste('Amel',c(1:6),sep ='_'))
files <- file.path('quants', aech_files, "quant.sf")
names(files) <- aech_files
tx2gene <- read.csv('Amel_gemoma_t2g.txt',header = F, sep = '\t')
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
batch = factor(c(rep(c('gyne','worker'),each =3 )))#,     c('minor','gyne','minor','gyne')))

sampleTable <- data.frame(condition = batch)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition)# + type )
dds <- DESeq(dds)
res <- results(dds,contrast = c('condition','gyne','worker'),alpha = 1e-3)
res <- results(dds,contrast = c('condition','gyne','worker'),alpha = 1e-3,lfcThreshold = log2(1.5))

res <- lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
#write.table(res, "Amel_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')

gyne_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange > 0),])
worker_bias = rownames(res[which(res$padj < 0.01 & res$log2FoldChange < 0),])
summary(res)
#ant_res = read.table('ant_res.txt',header = T,sep = '\t')
ant_res = read.table('../../ortholog/ortholog_gexpr_TPM_five_sp/candidate_deg.txt',header = T,sep = '\t')
t2gene = read.table('~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/ortholog/ortholog_gexpr_TPM_four_sp/candidate_gene/Aech_gemoma_t2g.txt')
library(readr)
tblast <- read_delim("~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/ortholog/ortholog_gexpr_TPM_four_sp/candidate_gene/proteins.fasta.blastout", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
t2gene$V1 = toupper(t2gene$V1)
t2gene$V2 = paste("Aech",t2gene$V2,sep = '_')
head(t2gene)
head(tblast)
t2gene$blast = tblast$X17[match(t2gene$V1,tblast$X1)]
t2gene$evalue = tblast$X15[match(t2gene$V1,tblast$X1)]
t2gene$coverage = tblast$X11[match(t2gene$V1,tblast$X1)]

gene_ortholog_table = read.table('../../ortholog/ortholog_ant_bee/gene_table.proteinortho', col.names = c('Aech','Sinv','Mpha','Lhum','Cbir','Dqua','Amel'),stringsAsFactors = F)
gene_ortholog_table_ant = read.table('../../ortholog/ortholog_seven/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua'),stringsAsFactors = F)
gene_ortholog_table  = gene_ortholog_table[which(gene_ortholog_table$Aech %in% gene_ortholog_table_ant$Aech),]
#gene_ortholog_table$Amel = paste('Amel',gene_ortholog_table$Amel,sep = '_')
res.data = data.frame(res)
res.data$Aech_ID = as.character(gene_ortholog_table$Aech)[match(row.names(res.data),gene_ortholog_table$Amel)]
res.data$annotation = t2gene$blast[match(res.data$Aech_ID,t2gene$V2)]
#res = data.frame(res)
table(res.data$padj < 1e-3 & !is.na(res.data$Aech_ID))
table(res.data$Aech_ID  %in%  rownames(ant_res) & res.data$padj < 1e-3 )

ant_res$log2_bee = res.data$log2FoldChange[match(rownames(ant_res),res.data$Aech_ID)]
ant_res$padj_bee = res.data$padj[match(rownames(ant_res),res.data$Aech_ID)]
tmp = ant_res[,c(3,1,4,2,5)]
names(tmp) = c("Gene_name",'Log2_ant','Log2_bee','padj_ant','padj_bee')
tmp = tmp[order(tmp$Log2_ant*tmp$Log2_bee,decreasing = T),]
tmp
(ant_res[which(ant_res$padj < 1e-2 &ant_res$padj_bee < 1e-3),])
ant_res[which(ant_res$padj < 1e-2 &ant_res$padj_bee < 1e-3 & ant_res$log2FoldChange*ant_res$log2_bee > 0),]
#write.table(tmp,'ant_bee_deg.tsv',quote = F,sep = '\t',row.names = F, col.names = T)
fisher.test(matrix(c(15,39-15,893-15,6964-893-39+15),nrow = 2))
