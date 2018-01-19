log2_deg_f = log2_deg[!apply(log2_deg,1,anyNA),]
log2_deg_f$KO = GO_annotation$KO[match(row.names(log2_deg_f),GO_annotation$V2)]
log2_deg_f = log2_deg_f[!is.na(log2_deg_f$KO),]
log2_deg_kegg = log2_deg_f$KO
log2_deg_f = as.matrix(log2_deg_f[!is.na(log2_deg_f$KO),c(1:7)])
log2_deg_cor = abs(cor(t((log2_deg_f[,c(1:7)])),method = 's'))

insulin = which(log2_deg_kegg %in% kegg.gs[122][[1]])
gyne_bs = which(rowMins(log2_deg_f[,c(1:7)]) > 0)
worker_bs = which(rowMaxs(log2_deg_f[,c(1:7)]) < 0)

ig = c(insulin, gyne_bs)
iw = c(insulin, worker_bs,gyne_bs)


heatmap(log2_deg_cor[worker_bs,worker_bs])
sum(log2_deg_cor[gyne_bs,gyne_bs])*2/(length(gyne_bs)*(length(gyne_bs)-1))
sum(log2_deg_cor[worker_bs,worker_bs])*2/(length(worker_bs)*(length(worker_bs)-1))
sample_row = sample(row.names(log2_deg_cor),length(worker_bs))
sum(log2_deg_cor[sample_row,sample_row])*2/(length(sample_row)*(length(sample_row)-1))

sum(log2_deg_cor[insulin,insulin])*2/(length(insulin)*(length(insulin)-1))
sum(log2_deg_cor[ig,ig])*2/(length(ig)*(length(ig)-1))
sum(log2_deg_cor[iw,iw])*2/(length(iw)*(length(iw)-1))

path_connet = data.frame(names = names(kegg.gs), connect = NA)
row.names(path_connet) = names(kegg.gs)
for(pathway in names(kegg.gs)){
  gene_set = kegg.gs[pathway][[1]]
  genes = which(log2_deg_kegg %in% gene_set)
  path_connet[pathway,2] = sum(log2_deg_cor[genes,genes])*2/(length(genes)*(length(genes)-1))
}






#####
log2_deg_cor_cum = apply(log2_deg_cor,1,FUN = function(x){
  sum(sort(abs(x),decreasing = T)[1:100])})

library(pathview)
library(gage)
kg.dme = kegg.gsets(species = 'ko',id.type = 'entrez')
kegg.gs=kg.dme$kg.sets#[kg.dme$sigmet.idx]
test = as.matrix(log2_deg_cor_cum)
row.names(test)  = GO_annotation$KO[match(names(log2_deg_cor_cum),GO_annotation$V2)]
test = test[!is.na(row.names(test)),]
tmp = gage(exprs = test,gsets = kegg.gs,same.dir=T,rank.test = T) #SAME.DIR can be used to test change with combined direction (kegg only)
head(tmp$greater,10)
tmp.sig = sigGeneSet(tmp,outname="tmp.sig")
require("RColorBrewer")
require("gplots")
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
col_breaks = c(seq(-5,-1,length=100),  # for red
               seq(-0.99,0.99,length=100),           # for yellow
               seq(1,5,length=100)) 
heatmap.2(tmp.sig$stats[,-c(1)],density.info="none",  trace="none",col=my_palette,dendrogram="none", Colv="NA",breaks=col_breaks,margins =c(10,20))
heatmap(-log10(rbind(tmp.sig$greater[,-c(1:5)])),scale = 'none') #Plotting the -log10(p value)
heatmap(tmp.sig$stats[,-c(1)],scale = 'none') #Plotting the stat
tmp.sig.up <- esset.grp(tmp$greater,exprs = test, gsets = kegg.gs, test4up = T, output = T, outname = "tmp.up", make.plot = T)
#####
go.dme=go.gsets(species="Fly",id.type = 'entrez')
go.bp=go.dme$go.sets[go.dme$go.subs$BP]
go.mf=go.dme$go.sets[go.dme$go.subs$MF]
go.cc=go.dme$go.sets[go.dme$go.subs$CC]
tmp = gage(exprs = test,gsets = go.bp,rank.test = T,FDR.adj = F,)
head(tmp$greater,n = 20)
head(tmp$greater[grep('light',rownames(tmp$greater)),])
tmp.sig = sigGeneSet(tmp,outname="tmp.sig")
tmp.sig
