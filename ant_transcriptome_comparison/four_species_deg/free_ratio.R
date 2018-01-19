w = read.table('free_ratio_clean.txt')
names(w) = c('geneID','Aech','Tcor','Ccos','Sinv','Mpha','Pbar','Lhum','Cbir','Dqua','Hsal')
w$median = rowMedians(as.matrix(w[,c(2:11)]))
log2_deg$var = apply(log2_deg[,c(2:5)],1,sd)
log2_deg$w = w$median[match(rownames(log2_deg),w$geneID)]
log2_deg_test = log2_deg[which(rownames(log2_deg) %in% rownames(res_sig)),]
plot(log2_deg$var, log2_deg$w)
cor.test(log2_deg_test$w,log2_deg_test$var,use = 'c',method = 's')

w$type = "NDE"
w$droID = t2gene$blast[match(w$geneID,t2gene$V2)]

w$type[w$geneID %in% rownames(deg_overlap_g)] = 'Gyne_bias'
w$type[w$geneID %in% rownames(deg_overlap_w)] = 'Worker_bias'

w$type[w$geneID %in% rownames(res_sig)] = 'Gyne_bias'
w$type[w$geneID %in% rownames(res_sig[which(res_sig$log2FoldChange< 0),])] = 'Worker_bias'
w$type[w$geneID %in% rownames(deg_overlap_nde)] = 'nde'
w$gene_dro = paste(w$geneID,w$droID,sep = '\n')
w$Mpha_log2 = log2_deg$Mpha_log2FoldChange[match(w$geneID,rownames(log2_deg))]

w$log2 = res$log2FoldChange[match(w$geneID,rownames(res))]
par(mfrow = c(1,1))
library(tidyr)
w_long <- gather(w, species, w, Aech:Hsal, factor_key=TRUE)
head(w_long)
colour_sp = c('blue','red')
library(ggplot2)
ggplot(data = w_long[which(w_long$type %in% c('Gyne_bias','Worker_bias')),],aes(x =gene_dro  , y  = w,  colour = type,label = species))+
  #geom_point(height = 0,width = 0.1)+
  #geom_hline(yintercept = median(w_long$w[w_long$type =='NDE']))+
  geom_hline(yintercept = colMedians(as.matrix(w[,c(2:11)])))+
#  geom_hline(yintercept = median(deg_overlap_nde$w,na.rm = T))+
  theme(axis.text.x=element_text(angle=70,hjust=1)) +
  ylim(c(0,2))+
  geom_text(size = 2)
 # geom_point()
 # scale_color_grey()+
  #scale_color_manual(values = colour_sp[c(1,1,1,1,1,1,1,0,0,0)+1])+
  ylim(c(0,2))
dim(w)
w_long[w_long$geneID == 'Aech_gene_11200',]
head(GO_annotation)
