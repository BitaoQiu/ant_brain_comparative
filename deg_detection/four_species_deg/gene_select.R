#Get the geneID for M.pharaonis and check the cDNA.
rownames(gene_ortholog_table) = gene_ortholog_table$Aech 
tmp2 = gene_ortholog_table[c(rownames(deg_overlap_w),rownames(deg_overlap_g)) ,c(1,3)]
head(tmp2)
tmp2
tmp4 = matrix(ncol = 10,nrow = dim(tmp2)[1])
rownames(tmp4) = tmp2$Aech
for(i in tmp2$Aech){tmp4[i,] = plotCounts(dds,normalized = T, gene=i, intgroup=c("caste",'species'), 
           returnData=TRUE)[c(17:26),1]}
tmp3 = rbind(deg_overlap_w,deg_overlap_g)[,c(2:5,6)]
head(tmp3)
ortholog_exp
tmp = ortholog_exp$abundance[rownames(ortholog_exp$abundance) %in% rownames(deg_overlap_g),c(seq(17,26,2),seq(18,26,2))]
tmp
head(GO_annotation)

assignment = read.csv("Mpha/assignment.tabular",header = T,sep = '\t')
head(assignment)
library(seqinr)
Mpha_fasta = read.fasta('Mpha/cds.fasta',as.string = T)
Mpha_fasta.data = data.frame(T_ID = names(Mpha_fasta),Seq = paste(Mpha_fasta))
output_tID = as.character(assignment[as.character(assignment$X.geneID) %in% gsub('Mpha_','',as.character(tmp2$Mpha)),'transcript'])
Mpha_out = Mpha_fasta.data[as.character(Mpha_fasta.data$T_ID) %in% output_tID,]
Mpha_out$gene_ID = assignment$X.geneID[match(Mpha_out$T_ID, as.character(assignment$transcript))]
Mpha_out$Seq = toupper(Mpha_out$Seq)
tmp_out = Mpha_out[,c(3,2)]
library(seqRFLP)
seq_out = dataframe2fas(tmp_out)
head(seq_out)
write.fasta(seq_out,file = 'Mpha/Mpha.fasta')
