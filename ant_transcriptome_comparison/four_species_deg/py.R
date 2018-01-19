library('RUVSeq')
library(vegan)
log2_deg_f = scale(log2_deg[!apply(log2_deg,1,anyNA),])
log2_deg_dist = dist(t(combat_edata))
#log2_deg_dist =  as.dist(1-cor(log2_deg_f,method = 's'))
par(mfrow = c(1,1))
plot(hclust(dist(t(log2_deg_f)), "average"))
plot(hclust(as.dist(1-cor(log2_deg_f,method = 's')), "average"))


dist_tmp = apply(combat_edata,1,FUN = function(x){
  mantel(log2_deg_dist,dist(x),method = 's',parallel = 4)$statistic})
  
sort(dist_tmp,decreasing = T)[1:10]
plot(hclust(dist(combat_edata['Aech_gene_11917',])))
plot(hclust(dist(log2_deg_f['Aech_gene_13852',])))
plot(hclust(dist(t(combat_edata))))
plot(hclust(as.dist(1-cor(combat_edata,method = 's'))))

#####
library(pathview)
library(gage)
kg.dme = kegg.gsets(species = 'dme',id.type = 'entrez')
kegg.gs=kg.dme$kg.sets[kg.dme$sigmet.idx]
test = as.matrix(dist_tmp)
row.names(test)  = GO_annotation$entrez[match(names(dist_tmp),GO_annotation$V2)]
test = test[!is.na(row.names(test)),]
tmp = gage(exprs = test,gsets = kegg.gs,same.dir=T,rank.test = T) #SAME.DIR can be used to test change with combined direction (kegg only)
tmp_greater = tmp$greater
plot(tmp_greater[,3],tmp_greater[,5],log = 'xy',xlab = 'p.value',ylab = 'Number of genes in pathway',pch = 4)
tmp_select = c(1:10)
text(tmp_greater[tmp_select,3], tmp_greater[tmp_select,5],
     labels=sub(pattern = 'dme......','',rownames(tmp_greater))[tmp_select], cex= .5,pos=4,col = 'red')
abline(v= 0.05,lty = 2)
head(tmp_greater)
#Pathview:
plot_path = kg.dme$sigmet.idx[grep('ko00020',names(kg.dme$kg.sets))]
pv.out <- pathview(gene.data = test, pathway.id = '00020', species = "ko", out.suffix = "test", kegg.native = T)
#####
go.dme=go.gsets(species="Fly",id.type = 'entrez')
go.bp=go.dme$go.sets[go.dme$go.subs$BP]
go.mf=go.dme$go.sets[go.dme$go.subs$MF]
go.cc=go.dme$go.sets[go.dme$go.subs$CC]
test = as.matrix(dist_tmp)
row.names(test)  = GO_annotation$entrez[match(names(dist_tmp),GO_annotation$V2)]
test = test[!is.na(row.names(test)),]
tmp = gage(exprs = test,gsets = go.bp,rank.test = T)

tmp_greater = tmp$greater
plot(tmp_greater[,3],tmp_greater[,5],log = 'xy',xlab = 'p.value',ylab = 'Number of genes in pathway',pch = 4)
tmp_select = c(1:10)
text(tmp_greater[tmp_select,3], tmp_greater[tmp_select,5],
     labels=sub(pattern = 'GO.........','',rownames(tmp_greater))[tmp_select], cex= .8,pos=4,col = 'red')
abline(v= 0.05,lty = 2)
head(tmp_greater[,c(3,5)], n  = 20)
#####
head(tmp$greater,n = 20)
#Tried my best, but no interesting signal? or 
#####
