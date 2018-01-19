#Had better do some filtering to make the ratio more robust.
library("pheatmap")
#filter = rownames(dds[which(rowSums(counts(dds)) > 8),])
filter = rownames(dds[which(rowMax(counts(dds)) > 5),])
abundance = txi.salmon$abundance[which(row.names(txi.salmon$abundance) %in% filter),] + 0.01
gene_age = read.table('../../synteny/Cbir_age.txt',header = T)
gene_age = read.table('../../synteny_seven//Cbir_age.txt',header = T)

head(tx2gene)
mpha_exp = abundance#assay(vsd)
#mpha_exp = scale(mpha_exp)
#mpha_exp = cbind(mpha_exp[,c(1,2)]- rowMeans(mpha_exp[,c(1,2)]),mpha_exp[,c(3,4)]- rowMeans(mpha_exp[,c(3,4)]),mpha_exp[,c(5,6)]- rowMeans(mpha_exp[,c(5,6)]),
#                 mpha_exp[,c(7,8)]- rowMeans(mpha_exp[,c(7,8)]))#,mpha_exp[,c(9,10)]- rowMeans(mpha_exp[,c(9,10)]))
#mpha_exp = mpha_exp[,sample(colnames(mpha_exp),8)]
#mpha_exp = log(mpha_exp)
mpha_exp = log2(mpha_exp[,c(1:4)]/abundance[,c(5:8)])#assay(vsd)
#mpha_exp = mpha_exp[!apply(mpha_exp,1,sd) == 0,]
#mpha_net = cor(t(mpha_exp),method= 's')
#####
#mpha_net_ac = apply(mpha_net,1,FUN = function(x){sum(sort(abs(x),decreasing = T)[1:1000])})
gene_age$geneID = tx2gene$V2[match(rownames(gene_age),tx2gene$V1)]

#gene_age$ac = mpha_net_ac[match(gene_age$geneID,names(mpha_net_ac))]/1000
gene_age = gene_age[which(gene_age$geneID %in% filter),]
#gene_age = gene_age[!is.na(gene_age$ac),]
gene_age$age = factor(gene_age$age,
                      levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
#plot(gene_age$age,gene_age$ac,pch = 20,ylim = c(0,1.1),xlab = 'Origin of genes',ylab = "Average connectivities")

#wilcox.test(gene_age$ac[gene_age$age == 'neoptera'],gene_age$ac[gene_age$age == 'TRG'])
table(gene_age$age)
######

#####
gene_age$caste = 'non-bias'
gene_age$caste[gene_age$geneID %in% gyne_bias] = 'gyne'
gene_age$caste[gene_age$geneID %in% worker_bias] = 'worker'
gene_age$caste = factor(gene_age$caste, levels = c('gyne','worker','non-bias'))
write.table(as.data.frame(table(gene_age$caste,gene_age$age)),'Cbir_age_odd.csv',sep = '\t',quote = F)
#boxplot(ac ~ caste, data = gene_age,ylim = c(0,1.1),xlab = 'Caste',ylab = "Average connectivities")
#plot(table( gene_age$age,gene_age$caste))
#####
#####
library("pheatmap")
library("RColorBrewer")
levels(gene_age$age) = c('Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Species specific')
levels(gene_age$caste) = c("Gyne bias","Worker bias","NDE")
n = 500
tmp_caste = c(as.character(sample(gene_age[which(gene_age$caste %in% c('Gyne bias','Worker bias')),'geneID'], n)))
#tmp_caste = c(as.character(sample(gene_age[which(gene_age$caste %in% c('Gyne bias')),'geneID'], n/2)),
#              as.character(sample(gene_age[which(gene_age$caste %in% c('Worker bias')),'geneID'], n/2)))
tmp = tmp_caste
age_bias = data.frame(Evolutionary.origin = gene_age$age[match(tmp, gene_age$geneID)],
                      Average.connectivity = apply(abs(mpha_net[tmp,tmp]),1,mean),
                      row.names = tmp)
caste_bias = data.frame("Caste.expression" = gene_age$caste[match(tmp_caste, gene_age$geneID)],
                        row.names = tmp_caste)
caste_bias = data.frame("Caste.expression" = res$log2FoldChange[match(tmp_caste, rownames(res))],
                        row.names = tmp_caste)
caste_bias$Caste.expression[abs(caste_bias$Caste.expression) >
                              max(caste_bias$Caste.expression)] = -max(caste_bias$Caste.expression)
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

pheatmap((mpha_net[tmp_caste,tmp_caste]), cluster_rows=T, show_rownames=F,show_colnames = F, cluster_cols=T,annotation_col=age_bias,annotation_row = caste_bias,
         annotation_colors = ann_colors,main = "L.humile",clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation')
#####
boxplot(ac ~ age, data = gene_age[gene_age$caste!= "NDE",])
abline(h = median(gene_age$ac[gene_age$age == 'Neoptera' & gene_age$caste!= "NDE"]))
boxplot(ac ~ age, data = gene_age)
abline(h = median(gene_age$ac[gene_age$age == 'Neoptera' ]))
caste_gene = as.character(gene_age[which(gene_age$caste %in% c('Gyne bias','Worker bias','NDE')),'geneID'])
caste_net = mpha_net[caste_gene,caste_gene]
caste_ac = data.frame(ac = apply(abs(caste_net),1,mean),
                      age = gene_age$age[match(caste_gene, gene_age$geneID)],
                      caste = droplevels(gene_age$caste[match(caste_gene, gene_age$geneID)]),
                      log2 = res$log2FoldChange[match(caste_gene, rownames(res))])
head(caste_ac)
boxplot(ac ~ age , data = caste_ac)
#write.table(caste_ac,file = 'Lhum_caste_ac.txt',sep = '\t',quote = F)
abline(h = median(caste_ac$ac[caste_ac$age == 'Neoptera']))
summary(lm(ac ~ age, data = caste_ac))
summary(lm(ac ~ log2, data = caste_ac))
summary(lm(ac ~ caste*age, data = caste_ac))
summary(lm(ac ~ caste, data = caste_ac))
cor.test(caste_ac$log2,caste_ac$ac,method = 'k')
#####
caste_gene_w = as.character(gene_age[which(gene_age$caste %in% c('Worker bias')),'geneID'])
caste_gene_g = as.character(gene_age[which(gene_age$caste %in% c('Gyne bias')),'geneID'])
caste_gene_b = as.character(gene_age[which(gene_age$caste %in% c('NDE')),'geneID'])
caste_net_w = abs(mpha_net[caste_gene_w,caste_gene_w])
mean(caste_net_w)
caste_net_g = abs(mpha_net[caste_gene_g,caste_gene_g])
mean(caste_net_g)
caste_net_b = abs(mpha_net[caste_gene_b,caste_gene_b])
mean(caste_net_b)
####
source("http://www.yilab.gatech.edu/pcor.R")
pcor.test(x=caste_ac$log2,y = caste_ac$ac,z =as.numeric(caste_ac$age), method = 's')
cor.test(caste_ac$ac,as.numeric(caste_ac$age),method = 'k')
cor.test(caste_ac$log2,as.numeric(caste_ac$age),method = 'k')
