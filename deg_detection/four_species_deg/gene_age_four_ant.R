#gene_age
gene_age = read.table('../../synteny/Aech_age.txt',header = T)
t2gene = read.table('../../deg_salmon_gemoma/Aech/Aech_gemoma_t2g.txt',header = F)
t2gene$V1 = toupper(t2gene$V1)
t2gene$V2 = paste('Aech_', t2gene$V2, sep = '')

gene_age$geneID = t2gene$V2[match(rownames(gene_age),t2gene$V1)]

tmp = data.frame(res[,c("log2FoldChange","padj")])

tmp = data.frame(log2_deg)

tmp$age = gene_age$age[match(rownames(tmp),gene_age$geneID)]
tmp$age[tmp$age %in% c("TRG",'myrmicinae')] = 'formicidae'
tmp$age = factor(tmp$age,levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae'))
levels(tmp$age) = c('Pre-neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae')
tmp$caste_bias = "NDE"
tmp$caste_bias[tmp$log2FoldChange > 0 &tmp$padj < 0.05] = 'Gyne'
tmp$caste_bias[tmp$log2FoldChange < 0 &tmp$padj < 0.05] = 'Worker'

tmp$caste_bias = "NDE"
tmp$caste_bias[which(rowSums(padj_deg[,c(2:5)] < .05) >= 4 & rowSums(log2_deg[,c(2:5)] > 0) == 4)] = 'Gyne'
tmp$caste_bias[which(rowSums(padj_deg[,c(2:5)] < .05) >= 4 & rowSums(log2_deg[,c(2:5)] < 0) == 4)] = 'Worker'

tmp$caste_bias = "NDE"
tmp[row.names(res[which(res$log2FoldChange>0 & res$padj < 0.01),]),'caste_bias'] = 'Gyne'
tmp[row.names(res[which(res$log2FoldChange<0 & res$padj < 0.01),]),'caste_bias'] = 'Worker'


table(tmp$caste_bias,tmp$age)
round(table(tmp$caste_bias,tmp$age)/rep(table(tmp$age),each = 3),2)

sum(tmp$caste_bias %in% c("Gyne",'Worker'))
