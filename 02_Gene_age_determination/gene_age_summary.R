library("pheatmap")
#Using the output of DESeq2 and combine with gene age information from age_formating.R
#Note: Here is only using M.pharaonis as an example, we did the same for other species.

filter = rownames(dds[which(rowMeans(counts(dds)) > 0),]) #Filter genes without expression level (maybe pseudo-genes)
abundance = txi.salmon$abundance[which(row.names(txi.salmon$abundance) %in% filter),] + 0.01
gene_age = read.table('Mpha_age.txt',header = T,stringsAsFactors = F)

gene_age$geneID = tx2gene$V2[match(rownames(gene_age),tx2gene$V1)]
gene_age = gene_age[which(gene_age$geneID %in% rownames(abundance)),]
gene_age$age = factor(gene_age$age,
                      levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
table(gene_age$age)
######
worker_bias = row.names(res[which(res$padj< 0.05 & res$log2FoldChange < 0),])
gyne_bias = row.names(res[which(res$padj< 0.05 & res$log2FoldChange > 0),])

gene_age$caste = 'non-bias'
gene_age$caste[gene_age$geneID %in% worker_bias] = 'worker'
gene_age$caste[gene_age$geneID %in% gyne_bias] = 'gyne'
gene_age$caste = factor(gene_age$caste, levels = c('gyne','worker','non-bias'))
write.table(as.data.frame(table(gene_age$caste,gene_age$age)),'Mpha_age_odd.csv',sep = '\t',quote = F)
