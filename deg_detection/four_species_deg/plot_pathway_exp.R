#log2_deg = log2_deg[!apply(log2_deg,1,anyNA),]
pathway_edata_c = pathway_edata[,seq(1,34,2)]-pathway_edata[,seq(2,34,2)]
pathway_edata_c = exp[,seq(1,34,2)]-exp[,seq(2,34,2)]
rownames(pathway_edata_c) = GO_annotation$entrez[match(rownames(pathway_edata_c),GO_annotation$V2)]
pathway_edata_c=pathway_edata_c[!apply(pathway_edata_c,1,anyNA),]
path_target = "dme04745 Phototransduction - fly"
tmp_path = pathway_edata_c[rownames(pathway_edata_c) %in% kegg.gs[path_target][[1]],]
tmp_path = pathway_edata_c[rownames(pathway_edata_c) %in% go.bp[path_target][[1]],]

tmp_path[tmp_path > abs(min(tmp_path,na.rm = T))] = abs(min(tmp_path,na.rm = T))
tmp_path[abs(tmp_path) > max(tmp_path,na.rm = T)] = -max(tmp_path,na.rm = T)
annotation_col = data.frame(row.names = colnames(pathway_edata_c), species = c(rep('A.echinatior',4),rep('L.humile',4),rep('M.pharaonis',5),rep('S.invicta',4)))
sp_color = grey.colors(4,start = 0.1,end = 1)
ann_colors = list(species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],L.humile = sp_color[3],M.pharaonis = sp_color[4]))
row_id = gsub('Dmel_','',GO_annotation$name[match(rownames(tmp_path),as.character(GO_annotation$entrez))])
n_break = 100
colors <- colorRampPalette(c(rgb(62,113,178, maxColorValue = 255),"white", rgb(215,50,40,maxColorValue = 255)),interpolate ='linear', alpha = TRUE)(n_break)
pheatmap(tmp_path,color = colors,scale = 'none',annotation_col = annotation_col,annotation_colors = ann_colors,
         show_colnames = F,labels_row = row_id,main = path_target)
quick_check(pathway_edata,sampleTable,path_target,gene_set = kegg.gs)
quick_check(pathway_edata,sampleTable,path_target,gene_set = go.bp)

