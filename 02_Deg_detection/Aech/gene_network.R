#Had better do some filtering to make the ratio more robust.
library("pheatmap")
filter = rownames(dds[which(rowMeans(counts(dds)) > 0),])
#filter = rownames(dds[which(rowMax(counts(dds)) > 5),])
abundance = txi.salmon$abundance[which(row.names(txi.salmon$abundance) %in% filter),] + 0.01
gene_age = read.table('../../synteny_seven/loss10//Aech_age.txt',header = T,stringsAsFactors = F)
gene_age = read.table('../../synteny_seven/loss20/Aech_age.txt',header = T,stringsAsFactors = F)
gene_age = read.table('../../synteny_seven/loss05//Aech_age.txt',header = T,stringsAsFactors = F)

gene_age$age = as.character(gene_age$age)
#blast_result = read.table('../blast_result/Aech.blast',header = F)
#blast_result$V1 = gsub('Aech_','', blast_result$V1)
#blast_result = blast_result[which(blast_result$V13 < 1e-5),]
#gene_age[which(rownames(gene_age) %in% blast_result$V1 & gene_age$age == 'TRG'),'age'] = 'homolog'
head(tx2gene)
library("pheatmap")
#mpha_exp = abundance#assay(vsd)
#mpha_exp = cbind(mpha_exp[,c(1:3)]- rowMeans(mpha_exp[,c(1:3)]),mpha_exp[,c(4:6)]- rowMeans(mpha_exp[,c(4:6)]),mpha_exp[,c(7:9)]- rowMeans(mpha_exp[,c(7:9)]),
#                 mpha_exp[,c(10:12)]- rowMeans(mpha_exp[,c(10:12)]))
#head(mpha_exp)
#mpha_exp = mpha_exp[!apply(mpha_exp,1,sd) == 0,]
#mpha_net = cor(t(mpha_exp),method= 's')
#####
gene_age$geneID = tx2gene$V2[match(rownames(gene_age),tx2gene$V1)]

gene_age = gene_age[which(gene_age$geneID %in% rownames(abundance)),]
gene_age$age = factor(gene_age$age,
                      levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
table(gene_age$age)
######
gene_age$caste = 'non-bias'
gene_age$caste[gene_age$geneID %in% minor_bias] = 'minor'
gene_age$caste[gene_age$geneID %in% major_bias] = 'major'
gene_age$caste[gene_age$geneID %in% gyne_bias] = 'gyne'
gene_age$caste[gene_age$geneID %in% worker_bias] = 'worker'

gene_age$caste = factor(gene_age$caste, levels = c('gyne','major','minor','worker','non-bias'))
write.table(as.data.frame(table(gene_age$caste,gene_age$age)),'../gene_age_five_10/Aech_age_odd.csv',sep = '\t',quote = F)
write.table(as.data.frame(table(gene_age$caste,gene_age$age)),'../gene_age_five_20/Aech_age_odd.csv',sep = '\t',quote = F)
write.table(as.data.frame(table(gene_age$caste,gene_age$age)),'../gene_age_five_05/Aech_age_odd.csv',sep = '\t',quote = F)
#####
library("pheatmap")
library("RColorBrewer")
levels(gene_age$age) = c('Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Species specific')
levels(gene_age$caste) = c("Gyne bias","Major bias",'Minor bias',"Worker bias","NDE")
n = 500
tmp_caste = as.character(sample(gene_age[which(gene_age$caste %in% c("Gyne bias","Major bias",'Minor bias',"Worker bias")),'geneID'], n))
tmp = tmp_caste
age_bias = data.frame(Evolutionary.origin = gene_age$age[match(tmp, gene_age$geneID)],
                      Average.connectivity = rowSums(abs(mpha_net[tmp,tmp]))/n,
                      row.names = tmp)
caste_bias = data.frame("Caste.expression" = gene_age$caste[match(tmp_caste, gene_age$geneID)],
                        row.names = tmp_caste)

tmp_col = gray.colors(8, start = 0.1, end = .9,gamma = 2.8)
ann_colors = list(
  Evolutionary.origin = c(Neoptera =tmp_col[1],Endopterygota = tmp_col[2],Hymenoptera=tmp_col[3],Aprocrita = tmp_col[4],Aculeata = tmp_col[5],
                          Formicidae = tmp_col[6],Myrmicinae = tmp_col[7],`Species specific` = tmp_col[8]),
  "Caste.expression" = c(`Gyne bias` = rgb(1,0,0,0.5),"Major bias" = rgb(0,1,0,0.5),'Minor bias'=rgb(1,0,1,0.5),`Worker bias` = rgb(0,0,1,0.5)),
  Average.connectivity = colorRampPalette( (brewer.pal(9, "Blues")) )(255))

pheatmap(abs(mpha_net[tmp_caste,tmp_caste]), cluster_rows=T, show_rownames=F,show_colnames = F, cluster_cols=T,annotation_col=age_bias,annotation_row = caste_bias,
         annotation_colors = ann_colors,main = "A.echinatior")
#####

#####
caste_gene = as.character(gene_age[which(gene_age$caste %in% c("Gyne bias","Major bias",'Minor bias',"Worker bias")),'geneID'])
caste_net = mpha_net[caste_gene,caste_gene]
caste_ac = data.frame(ac = apply(abs(caste_net),1,mean),
                      age = gene_age$age[match(caste_gene, gene_age$geneID)],
                      caste = droplevels(gene_age$caste[match(caste_gene, gene_age$geneID)]))
head(caste_ac)
boxplot(ac ~ age+caste, data = caste_ac)
abline(h = median(caste_ac$ac[caste_ac$age == 'Neoptera']))
summary(lm(ac ~ caste, data = caste_ac))
summary(lm(ac ~ age, data = caste_ac))
summary(lm(ac ~ age*caste, data = caste_ac))
