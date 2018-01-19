#PY 2, pathway based
library('RUVSeq')
library(vegan)
log2_deg_dist = as.dist(1-cor(combat_edata,method = 's'))
pathway_edata = combat_edata
rownames(pathway_edata) = GO_annotation$entrez[match(rownames(pathway_edata),GO_annotation$V2)]

pathway_p_value = data.frame(row.names = names(kegg.gs),N_gene = rep(NA,length(names(kegg.gs))), cor = rep(NA,length(names(kegg.gs))),
                             pvalue = rep(NA,length(names(kegg.gs))))
for(i in names(kegg.gs)){
  tmp_data = pathway_edata[rownames(pathway_edata) %in% kegg.gs[i][[1]],]
  if(!is.null(dim(tmp_data))){
    pathway_p_value[i,'N_gene'] = dim(tmp_data)[1]
    pathway_p_value[i,'cor'] = try(mantel(log2_deg_dist, as.dist(1-cor(tmp_data,method = 's')),method = 's',parallel = 4)$statistic,silent = T)
    pathway_p_value[i,'pvalue'] = try(mantel(log2_deg_dist, as.dist(1-cor(tmp_data,method = 's')),method = 's',parallel = 4)$signif,silent = T)
  }
}
pathway_p_value_filter = pathway_p_value[which((pathway_p_value$N_gene > 8 & pathway_p_value$cor > 0)),]
pathway_p_value_filter$cor = round(as.numeric(pathway_p_value_filter$cor),digits = 3)
path_tmp_table = pathway_p_value_filter[order(pathway_p_value_filter$cor,decreasing = T)[c(1:20)],]
path_tmp_table
quick_check(pathway_edata,sampleTable,"dme03008 Ribosome biogenesis in eukaryotes",gene_set = kegg.gs)

plot(pathway_p_value_filter$cor,pathway_p_value_filter$N_gene,log = 'xy',xlab = 'Correlation',ylab = 'Number of genes in pathway',pch = '.')
select = order(pathway_p_value_filter$cor,decreasing = T)[c(1:20)]
text(pathway_p_value_filter[select,2], pathway_p_value_filter[select,1],
     labels=sub(pattern = 'dme......','',rownames(pathway_p_value_filter))[select], cex= .5,pos=1,col = 'red')
pathway_edata[rownames(pathway_edata) %in% kegg.gs["dme04080 Neuroactive ligand-receptor interaction"][[1]],]
log2_deg$entrez = GO_annotation$entrez[match(rownames(log2_deg),GO_annotation$V2)]
log2_deg$name = GO_annotation$name[match(rownames(log2_deg),GO_annotation$V2)]

log2_deg[log2_deg$entrez %in% kegg.gs["dme04080 Neuroactive ligand-receptor interaction"][[1]],c(2:5,9)]
#####
go_p_value = data.frame(row.names = names(go.bp),N_gene = rep(NA,length(names(go.bp))), cor = rep(NA,length(names(go.bp))),
                             pvalue = rep(NA,length(names(go.bp))))
for(i in names(go.bp)){
  tmp_data = pathway_edata[rownames(pathway_edata) %in% go.bp[i][[1]],]
  if(!is.null(dim(tmp_data))){
    go_p_value[i,'N_gene'] = dim(tmp_data)[1]
    go_p_value[i,'cor'] = try(mantel(log2_deg_dist, as.dist(1-cor(tmp_data,method = 's')),method = 's',parallel = 4)$statistic,silent = T)
    go_p_value[i,'pvalue'] = try(mantel(log2_deg_dist, as.dist(1-cor(tmp_data,method = 's')),method = 's',parallel = 4)$signif,silent = T)
  }
}
go_p_value_filter = go_p_value[which((go_p_value$N_gene > 8 )),]
go_p_value_filter$cor = round(as.numeric(go_p_value_filter$cor ),digits = 3)
go_p_value_filter[order(go_p_value_filter$cor,decreasing = T)[c(1:20)],]
plot(go_p_value_filter$cor,go_p_value_filter$N_gene,log = 'xy',xlab = 'Correlation',ylab = 'Number of genes in pathway',pch = 4)
select = order(go_p_value_filter$cor,decreasing = T)[c(1:10)]
text(go_p_value_filter[select,2], go_p_value_filter[select,1],
     labels=sub(pattern = 'GO.........','',rownames(go_p_value_filter))[select], cex= .5,pos=2,col = 'red')


path_tmp = pathway_edata[rownames(pathway_edata) %in% go.bp["GO:0008306 associative learning"][[1]],]
path_tmp = pathway_edata[rownames(pathway_edata) %in% kegg.gs["dme01230 Biosynthesis of amino acids"][[1]],]

quick_check(pathway_edata,sampleTable,"dme04068 FoxO signaling pathway",gene_set = kegg.gs)
