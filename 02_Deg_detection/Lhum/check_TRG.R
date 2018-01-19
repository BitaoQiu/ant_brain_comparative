blast_result = read.table('Lhum.blast')
head(blast_result)
blast_result = blast_result[which(blast_result$V15 > 50),]

table(rownames(gene_age) %in% blast_result$V1,gene_age$age)
table(rownames(gene_age) %in% blast_result$V1,gene_age$caste)
table(gene_age$caste,gene_age$age)
