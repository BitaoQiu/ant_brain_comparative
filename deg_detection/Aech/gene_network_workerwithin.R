#Had better do some filtering to make the ratio more robust.
library("pheatmap")
filter = rownames(dds[which(rowSums(counts(dds)) > 12),])
abundance = txi.salmon$abundance[which(row.names(txi.salmon$abundance) %in% filter),] + 0.01
gene_age = read.table('../../synteny/Aech_age.txt',header = T)
head(tx2gene)
mpha_exp = abundance#assay(vsd)
mpha_exp = cbind(mpha_exp[,c(2,3)]- rowMeans(mpha_exp[,c(2,3)]),mpha_exp[,c(5,6)]- rowMeans(mpha_exp[,c(5,6)]),mpha_exp[,c(8,9)]- rowMeans(mpha_exp[,c(8,9)]),
                 mpha_exp[,c(11,12)]- rowMeans(mpha_exp[,c(11,12)]))
mpha_exp = log2(abundance[,seq(2,12,3)]/abundance[,seq(3,12,3)])#assay(vsd)
#mpha_exp = log2(abundance[,seq(1,12,3)]/abundance[,seq(2,12,3)])#assay(vsd)
mpha_exp = mpha_exp[!apply(mpha_exp,1,sd) == 0,]
mpha_net = cor(t(mpha_exp),method= 's')
#####
#mpha_net_ac = apply(mpha_net,1,FUN = function(x){sum(sort(abs(x),decreasing = T)[1:1000])})
gene_age$geneID = tx2gene$V2[match(rownames(gene_age),tx2gene$V1)]
gene_age = gene_age[which(gene_age$geneID %in% rownames(mpha_exp)),]
#gene_age$ac = mpha_net_ac[match(gene_age$geneID,names(mpha_net_ac))]/1000

#gene_age = gene_age[!is.na(gene_age$ac),]
gene_age$age = factor(gene_age$age,
                      levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
######
gene_age$caste = 'non-bias'
gene_age$caste[gene_age$geneID %in% major_bias_2] = 'major'
gene_age$caste[gene_age$geneID %in% minor_bias_2] = 'minor'

gene_age$caste = factor(gene_age$caste, levels = c('major','minor','non-bias'))
#####
#####
library("pheatmap")
library("RColorBrewer")
levels(gene_age$age) = c('Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Species specific')
levels(gene_age$caste) = c("Major bias","Minor bias","NDE")
n = 500
tmp_caste = as.character(sample(gene_age[which(gene_age$caste %in% c("Major bias","Minor bias","NDE")),'geneID'], n))
tmp = tmp_caste
age_bias = data.frame(Evolutionary.origin = gene_age$age[match(tmp, gene_age$geneID)],
                      Average.connectivity = rowSums(abs(mpha_net[tmp,tmp]))/n,
                      row.names = tmp)
caste_bias = data.frame("Caste.expression" = gene_age$caste[match(tmp_caste, gene_age$geneID)],
                        row.names = tmp_caste)
caste_bias = data.frame("Caste.expression" = res_bs$log2FoldChange[match(tmp_caste, rownames(res_bs))],
                        row.names = tmp_caste)
caste_bias$Caste.expression[caste_bias$Caste.expression >
                              abs(min(caste_bias$Caste.expression))] = abs(min(caste_bias$Caste.expression))
colnames(caste_bias) = "log2(caste expression ratio)"
tmp_col = gray.colors(8, start = 0.1, end = .9,gamma = 2.8)
ann_colors = list(
  Evolutionary.origin = c(Neoptera =tmp_col[1],Endopterygota = tmp_col[2],Hymenoptera=tmp_col[3],Aprocrita = tmp_col[4],Aculeata = tmp_col[5],
                          Formicidae = tmp_col[6],Myrmicinae = tmp_col[7],`Species specific` = tmp_col[8]),
  # "Caste.expression" = c(`Gyne bias` = rgb(1,0,0,0.5),`Worker bias` = rgb(0,0,1,0.5),NDE = rgb(0,0,0,.1)),
  "log2(caste expression ratio)" = colorRampPalette(c(rgb(62,113,178, maxColorValue = 255),"white", rgb(215,50,40,maxColorValue = 255)),interpolate ='linear', alpha = TRUE)(100),
  Average.connectivity = colorRampPalette( (brewer.pal(9, "Blues")) )(255))

pheatmap(abs(mpha_net[tmp_caste,tmp_caste]), cluster_rows=T, show_rownames=F,show_colnames = F, cluster_cols=T,annotation_col=age_bias,annotation_row = caste_bias,
         annotation_colors = ann_colors,main = "A.echinatior")
#####
#####
caste_gene = as.character(gene_age[which(gene_age$caste %in% c("Major bias","Minor bias","NDE")),'geneID'])
caste_net = mpha_net[caste_gene,caste_gene]
caste_ac = data.frame(ac = apply(abs(caste_net),1,mean),
                      age = gene_age$age[match(caste_gene, gene_age$geneID)],
                      caste = droplevels(gene_age$caste[match(caste_gene, gene_age$geneID)]),
                      log2 = res_bs$log2FoldChange[match(caste_gene, rownames(res_bs))])
head(caste_ac)
boxplot(ac ~ age, data = caste_ac)
abline(h = median(caste_ac$ac[caste_ac$age == 'Neoptera']))
summary(lm(ac ~ caste, data = caste_ac))
summary(lm(ac ~ log2, data = caste_ac))
summary(lm(ac ~ age, data = caste_ac))
summary(lm(ac ~ age*log2, data = caste_ac))
summary(lm(log2 ~ age, data = caste_ac))
boxplot(log2 ~ age, data = caste_ac)
abline(h = 0)
