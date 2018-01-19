library("RColorBrewer")
library("pheatmap")
deg_input = function(file_name, species, col_header){
  output_file = read.table(file_name)
  rownames(output_file) = paste(species, rownames(output_file), sep = '_')
  colnames(output_file) = paste(col_header, colnames(output_file), sep = '_')
  return(output_file)}

Aech_gyne_major_deg = deg_input('../../deg_ballgown/Aech_ERCC/Aech_gyne_major.csv','Aech','Aech_g_major')
Aech_gyne_minor_deg = deg_input('../../deg_ballgown/Aech_ERCC/Aech_gyne_minor.csv','Aech','Aech_g_minor')
Lhum_deg = deg_input('../../deg_ballgown/Lhum_ERCC//Lhum_gyne_worker.csv','Lhum','Lhum')
Mpha_deg = deg_input('../../deg_ballgown/Mpha_ERCC///Mpha_gyne_worker.csv','Mpha','Mpha')
Sinv_deg = deg_input('../../deg_ballgown/Sinv///Sinv_gyne_worker.csv','Sinv','Sinv')
Cbir_deg = deg_input('../../deg_ballgown/Cbir/Cbir_queen_worker.csv','Cbir','Cbir')
Dqua_deg = deg_input('../../deg_ballgown/Dqua/Dqua_queen_worker.csv','Dqua','Dqua')

ortholog_deg = cbind(Aech_gyne_major_deg[match(gene_ortholog_table$Aech, rownames(Aech_gyne_major_deg)),],
                     Aech_gyne_minor_deg[match(gene_ortholog_table$Aech, rownames(Aech_gyne_minor_deg)),],
                     Mpha_deg[match(gene_ortholog_table$Mpha, rownames(Mpha_deg)),],
                     Sinv_deg[match(gene_ortholog_table$Sinv, rownames(Sinv_deg)),],
                     Lhum_deg[match(gene_ortholog_table$Lhum, rownames(Lhum_deg)),],
                     Cbir_deg[match(gene_ortholog_table$Cbir, rownames(Cbir_deg)),],
                     Dqua_deg[match(gene_ortholog_table$Dqua, rownames(Dqua_deg)),])
log2_deg = ortholog_deg[,grep('log2Fold', colnames(ortholog_deg))]
padj_deg = ortholog_deg[,grep('padj', colnames(ortholog_deg))]
stat_deg = ortholog_deg[,grep('stat', colnames(ortholog_deg))]
log2_tmp = cbind(apply(log2_deg[,c(1,2)],1,min),log2_deg[,c(3:7)])
colnames(log2_tmp) = c('Aech','Mpha',"Sinv","Lhum","Cbir","Dqua")
library(limma)
a = vennCounts(padj_deg[,c(1:5)] < .1 & log2_deg[,c(1:5)] < 0)
vennDiagram(a)
log2_deg = ortholog_deg[,seq(4,42,6)]
plot(hclust(as.dist(1-cor(stat_deg, method = 's', use = 'complete.obs'))))

sampleDists = as.dist(1-cor(log2_tmp,method = 's', use = 'complete.obs'))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(log2_tmp)
colnames(sampleDistMatrix) <- colnames(log2_tmp)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,cluster_rows = T,cluster_cols = T,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#####

deg_union = rownames(padj_deg[which(rowSums(padj_deg[,c(2:5)] < .1,na.rm = T) > 0),])
#deg_union = rownames(log2_deg[which(rowSums(abs(log2_deg) > 1,na.rm = T) > 0),])

#####
gyne_deg = padj_deg < .1 & log2_deg > 0
head(gyne_deg)
gyne_deg_share = matrix(nrow = 7, ncol = 8,dimnames = list(colnames(gyne_deg),c(colnames(gyne_deg),'number')))
for(i in c(1:7)){
  for (j in c(1:7)){
    gyne_deg_share[i,j] =sum(gyne_deg[,i] & gyne_deg[,j],na.rm = T)/sum( gyne_deg[,i],na.rm = T)
  }
}
gyne_deg_share[,8] = colSums( gyne_deg,na.rm = T)
round(gyne_deg_share,3)

worker_deg = padj_deg < .1 & log2_deg < 0
head(worker_deg)
worker_deg_share = matrix(nrow = 7, ncol = 8,dimnames = list(colnames(worker_deg),c(colnames(gyne_deg),'number')))
for(i in c(1:7)){
  for (j in c(1:7)){
    worker_deg_share[i,j] =sum(worker_deg[,i] & worker_deg[,j],na.rm = T)/sum( worker_deg[,i],na.rm = T)
  }
}
worker_deg_share[,8] = colSums( worker_deg,na.rm = T)
round(worker_deg_share,3)

all_deg_share = matrix(nrow = 7, ncol = 8,dimnames = list(colnames(worker_deg),c(colnames(worker_deg),'Number')))
for(i in c(1:7)){
  for (j in c(1:7)){
    all_deg_share[i,j] =sum((worker_deg[,i] & worker_deg[,j]) |
                              (gyne_deg[,i] & gyne_deg[,j]),na.rm = T)/sum( worker_deg[,i]|gyne_deg[,i],na.rm = T)
  }
}
all_deg_share[,8] = colSums( worker_deg|gyne_deg,na.rm = T)
round(all_deg_share,3)


deg_union = rownames(padj_deg[which(rowSums(padj_deg < .1,na.rm = T) >1),])
length(deg_union)


######
#####
log2_deg = log2_deg[!apply(log2_deg,1,anyNA),]
tmp = t(log2_deg[sample(rownames(log2_deg),1000),])

sampleDists = as.dist(1-abs(cor(tmp,method = 's', use = 'everything')))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,#cluster_rows = F,cluster_cols = F,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
  