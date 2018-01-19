GO_annotation = read.table('../../ortholog/ortholog_gexpr_TPM/functional/Aech_annotation.tsv',header = T)
head(GO_annotation)
#####
#####
library(pathview)
library(gage)
#kg.dme = kegg.gsets(species = 'ko',id.type = 'entrez')
kg.dme = kegg.gsets(species = 'dme',id.type = 'entrez')

kegg.gs=kg.dme$kg.sets[kg.dme$sigmet.idx]
test = as.matrix(log2_deg)#It includes all genes, not just DEGs
row.names(test)  = GO_annotation$entrez[match(rownames(log2_deg),GO_annotation$V2)]
test = test[!is.na(row.names(test)),]
test_select = test[c(which(apply(test[,c(2:5)],1,min) > 0),which(apply(test[,c(2:5)],1,max) < 0)),]
colnames(test_select) = c('Aech_B','Aech_S','Mpha','Sinv','Lhum','Cbir','Dqua')
tmp = gage(exprs = pathway_edata_c,gsets = kegg.gs,same.dir=F,rank.test = T) #SAME.DIR can be used to test change with combined direction (kegg only)
head(tmp$greater,20)
tmp.sig = sigGeneSet(tmp,outname="tmp.sig")
tmp.sig
tmp.sig$greater[apply(tmp.sig$greater[,c(6:9)],1,max) < 0.1,]
tmp.sig$less[apply(tmp.sig$less[,c(6:9)],1,max) < 0.1,]

#annotation_col = data.frame(row.names = colnames(test), species = c('Aech','Aech','Mpha','Sinv','Lhum','Cbir','Dqua'),
#                            comparison = c("Gyne-Major","Gyne-Minor",rep("Gyne-Worker",5)))
annotation_col = data.frame(row.names = colnames(test), species = c('A.echinatior','A.echinatior','M.pharaonis','S.invicta','L.humile','Cbir','Dqua'))
select_ko = row.names(tmp.sig$stats[c(order(tmp.sig$stats[,1],decreasing = T)[1:6],order(tmp.sig$stats[,1],decreasing = F)[1:6]),])
select_ko = row.names(path_tmp_table)
sp_color = grey.colors(4,start = 0.1,end = 1)
ann_colors = list(species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],L.humile = sp_color[3],M.pharaonis = sp_color[4]))
n_break = 100
colors <- colorRampPalette(c(rgb(62,113,178, maxColorValue = 255),"white", rgb(215,50,40,maxColorValue = 255)),interpolate ='linear', alpha = TRUE)(n_break)

pheatmap(t(tmp$stats[which(rownames(tmp$stats) %in% select_ko)[1:12],-c(1)]),show_rownames =F,show_colnames = T,scale = 'none',annotation_row = annotation_col,
         annotation_colors = ann_colors,col=colors,breaks =seq(-10,10,length.out = n_break))
pheatmap(t(tmp.sig$stats[,-c(1)]),show_rownames =F,show_colnames = T,scale = 'none',annotation_row = annotation_col,
         annotation_colors = ann_colors,col=colors,breaks =seq(-10,10,length.out = n_break),rot = 90)
pheatmap(tmp.sig$stats[,-c(1)],show_rownames =T,show_colnames = F,scale = 'none',annotation_col = annotation_col,
         annotation_colors = ann_colors,col=colors,breaks =seq(-10,10,length.out = n_break) )
names(kg.dme$kg.sets)[grep('Insulin',names(kg.dme$kg.sets))][2]

plot_path = '04080'
pv.out <- pathview(gene.data = test[,c(2:5)], pathway.id = plot_path, species = "dme", kegg.native = T)
#####
test = as.matrix(log2_deg)
row.names(test)  = GO_annotation$entrez[match(rownames(log2_deg),GO_annotation$V2)]
test = test[!is.na(row.names(test)),]
colnames(test) = c('Aech_B','Aech_S','Mpha','Sinv','Lhum','Cbir','Dqua')
go.dme=go.gsets(species="Fly",id.type = 'entrez')
go.bp=go.dme$go.sets[go.dme$go.subs$BP]
go.mf=go.dme$go.sets[go.dme$go.subs$MF]
go.cc=go.dme$go.sets[go.dme$go.subs$CC]
tmp = gage(exprs = test[,c(2:5)],gsets = go.bp,same.dir = T,rank.test = T)
head(tmp$less,n = 20)
tmp.sig = sigGeneSet(tmp,outname="tmp.sig")
select_go = c(order(tmp.sig$stats[,1],decreasing = T)[1:10],order(tmp.sig$stats[,1],decreasing = F)[1:10])
select_go = which(abs(rowSums(tmp.sig$stats[,c(-1)]/abs(tmp.sig$stats[,c(-1)])) ) ==4)
select_go_order = c(order(tmp.sig$stats[select_go,1],decreasing = T)[1:10],order(tmp.sig$stats[select_go,1],decreasing = F)[1:10])

annotation_col = data.frame(row.names = colnames(test), species = c('A.echinatior','A.echinatior','M.pharaonis','S.invicta','L.humile','Cbir','Dqua'))
sp_color = grey.colors(4,start = 0.1,end = 1)
ann_colors = list(species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],L.humile = sp_color[3],M.pharaonis = sp_color[4]))
n_break = 100
colors <- colorRampPalette(c(rgb(62,113,178, maxColorValue = 255),"white", rgb(215,50,40,maxColorValue = 255)),interpolate ='linear', alpha = TRUE)(n_break)

pheatmap(tmp.sig$stats[select_go[select_go_order],-c(1)],show_rownames =T,show_colnames = F,scale = 'none',annotation_col = annotation_col,
         annotation_colors = ann_colors,col=colors,breaks =seq(-10,10,length.out = n_break) ,
         labels_row = sub(pattern = 'GO.........','',rownames(tmp.sig$stats[select_go[select_go_order],])))
pheatmap(tmp.sig$stats[,-c(1)],show_rownames =F,show_colnames = F,scale = 'none',annotation_col = annotation_col,
         annotation_colors = ann_colors,col=colors,breaks =seq(-10,10,length.out = n_break) ,
         labels_row = sub(pattern = 'GO.........','',rownames(tmp.sig$stats)))


tmp = gage(exprs = pathway_edata_c,gsets = go.bp,same.dir = T,rank.test = F)
head(tmp$less,n = 20)
tmp.sig = sigGeneSet(tmp,outname="tmp.sig")
