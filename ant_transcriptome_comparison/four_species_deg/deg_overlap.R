head(log2_deg)
log2_deg = ortholog_deg[,grep('log2Fold', colnames(ortholog_deg))]
padj_deg = ortholog_deg[,grep('padj', colnames(ortholog_deg))]
stat_deg = ortholog_deg[,grep('stat', colnames(ortholog_deg))]
basemean_deg = ortholog_deg[,grep('baseMean', colnames(ortholog_deg))]
library(limma)
colnames(padj_deg)[c(2:5)] = c('A.echinatior','M.pharaonis',"S.invicta",'L.humile')
a = vennCounts(padj_deg[,c(2:5)] < .05 & log2_deg[,c(2:5)] < 0)
a = vennCounts(padj_deg[,c(2:5)] < .05 & log2_deg[,c(2:5)] > 0)
vennDiagram(object = a,circle.col = c('red','lightgreen','purple','darkgreen'),main = 'Genes with worker-biased expression')
vennDiagram(object = a,circle.col = c('red','lightgreen','purple','darkgreen'),main = 'Genes with gyne-biased expression')

dim(padj_deg)
deg_overlap_w = round(padj_deg[which(rowSums(padj_deg[,c(2:5)] < .1) >= 4 & rowSums(log2_deg[,c(2:5)] < -.2) == 4),c(1:5)],4)
deg_overlap_w$dmel = GO_annotation$name[match(rownames(deg_overlap_w),GO_annotation$V2)]
deg_overlap_w
deg_overlap_g = round(padj_deg[which(rowSums(padj_deg[,c(2:5)] < .1) >= 4 & rowSums(log2_deg[,c(2:5)] > .2) == 4),c(1:5)],4)
deg_overlap_g$dmel = GO_annotation$name[match(rownames(deg_overlap_g),GO_annotation$V2)]
deg_overlap_g[,c(-1)]
deg_overlap_nde = round(padj_deg[which(rowSums(padj_deg[,c(1:7)] < .1) == 0),c(1:5)],4)
basemean_deg[rownames(deg_overlap_g),c(2:5)]
basemean_deg[rownames(deg_overlap_w),c(2:5)]

check_exp = function(ortholog_exp, gene_id,title =gene_id ){
  target_gene_exp = ortholog_exp[gene_id,]
  target_gene = data.frame(gyne = target_gene_exp[seq(1,34,2)],worker = target_gene_exp[seq(2,34,2)], 
                           species = c(rep(c('A.echinatior','L.humile'),each = 4),rep("M.pharaonis",5),rep("S.invicta",4)))
  ggplot(data =target_gene, aes(x = worker, y = gyne, colour = species) )+
    geom_point()+
    geom_abline(slope = 1, intercept = 0,lty = 2)+
    scale_x_sqrt(limits = c(0,max(target_gene_exp)+1))+
    scale_y_sqrt(limits = c(0,max(target_gene_exp)+1))+
    ggtitle(title)+xlab("Worker Expression (TPM)")+ylab("Gyne Expression (TPM)")}
t2gene[grep('insulin',t2gene$blast),]
check_exp(ortholog_exp,'Aech_gene_15120',title = 'GABA receptor')

check_exp_list = function(ortholog_exp, gene_id_list, title, gene_name=gene_id_list){
  target_gene_exp = ortholog_exp[gene_id_list,]
  target_gene = data.frame(gene = rep(gene_name, each = 17), 
                           gyne = c(t(target_gene_exp[,seq(1,34,2)])),
                           worker = c(t(target_gene_exp[,seq(2,34,2)])),
                           species = rep(c(rep(c('A.echinatior','L.humile'),each = 4),rep("M.pharaonis",5),rep("S.invicta",4)),length(gene_id_list)))
  ggplot(data =target_gene, aes(x = worker, y = gyne, colour = species) )+
    geom_point()+
    geom_abline(slope = 1, intercept = 0,lty = 2)+
#    scale_x_sqrt(limits = c(0,max(target_gene_exp)+1))+
 #   scale_y_sqrt(limits = c(0,max(target_gene_exp)+1))+
    facet_wrap(~gene,ncol = 2,scales = 'free')+
    theme(legend.position = 'bottom')+
    ggtitle(title)+xlab("Worker Expression (TPM)")+ylab("Gyne Expression (TPM)")}
check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_9986','Aech_gene_7060'),title = 's')
