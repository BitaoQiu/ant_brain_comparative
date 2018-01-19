pathway_edata = combat_edata
rownames(pathway_edata) = GO_annotation$entrez[match(rownames(pathway_edata),GO_annotation$V2)]

pathway_info = data.frame(row.names = names(kegg.gs),
                          N_origin = rep(NA,length(names(kegg.gs))), 
                          N_gene = rep(NA,length(names(kegg.gs))), 
                          bias_value = rep(NA,length(names(kegg.gs))),
                          sd = rep(NA,length(names(kegg.gs))),
                          cut_off = rep(NA,length(names(kegg.gs))),
                          cut_off_N = rep(NA,length(names(kegg.gs))),
                          abs_bias_value = rep(NA,length(names(kegg.gs))),
                          pvalue = rep(NA,length(names(kegg.gs))),
                          abs_pvalue = rep(NA,length(names(kegg.gs))))

#tmp_data = pathway_edata[rownames(pathway_edata) %in% kegg.gs[1][[1]],]
for(i in names(kegg.gs)){
  tmp_data = pathway_edata[rownames(pathway_edata) %in% kegg.gs[i][[1]],]
  if(!is.null(dim(tmp_data))){
    if (dim(tmp_data)[1] > 1){
      tmp_data_2 = tmp_data[,seq(1,34,2)] - tmp_data[,seq(2,34,2)]
      pathway_info[i,'N_origin'] = length(kegg.gs[i][[1]])
      pathway_info[i,'N_gene'] = dim(tmp_data_2)[1]
      pathway_info[i,'bias_value'] = sum(rowSums(tmp_data_2))/dim(tmp_data_2)[1]
      pathway_info[i,'sd'] = sum(rowSds(tmp_data_2))/dim(tmp_data_2)[1]
      pathway_info[i,'cut_off'] =  (sum(apply(tmp_data_2,1,min) > 0) + sum(apply(tmp_data_2,1,max) < 0))/dim(tmp_data_2)[1]
      pathway_info[i,'cut_off_N'] =  (sum(apply(tmp_data_2,1,min) > 0) + sum(apply(tmp_data_2,1,max) < 0))
      pathway_info[i,'abs_bias_value'] = sum(abs(rowSums(tmp_data_2)))/dim(tmp_data_2)[1]
      pathway_info[i,'pvalue'] = t.test(rowSums(tmp_data_2)/dim(tmp_data_2)[2])$p.value
      pathway_info[i,'abs_pvalue'] = t.test(abs(rowSums(tmp_data_2))/dim(tmp_data_2)[2])$p.value}
  }
}
pathway_info$pvalue = p.adjust(pathway_info$pvalue,method = 'fdr')
pathway_info$abs_pvalue = p.adjust(pathway_info$abs_pvalue,method = 'fdr')
#pathway_info = pathway_info[pathway_info$N_gene > 10,]
pathway_info[order(pathway_info$bias_value,decreasing = T)[1:20],]
pathway_info[order(pathway_info$bias_value,decreasing = F)[1:20],]
out = pathway_info[order(pathway_info$abs_bias_value,decreasing = T)[1:30],]
pathway_info[order(pathway_info$cut_off,decreasing = T)[1:20],]


pathway_info[order(pathway_info$pvalue,decreasing = F)[1:20],]
pathway_info[order(pathway_info$abs_pvalue,decreasing = F)[1:10],]
