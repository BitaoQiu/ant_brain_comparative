# Continued with the data of without_normalization.R
library(sva)
ortholog_exp_filtered = ortholog_exp[!apply(ortholog_exp,1, anyNA),] 

filter_table = sampleTable
batch = droplevels(filter_table$species)
# batch = droplevels(filter_table$colony) #Remove the # when normalising for colony identity.
modcombat = model.matrix(~1, data=filter_table)
combat_edata = ComBat(dat=exp, batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)

sampleDists = cor(combat_edata,method = 's')
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(combat_edata)
colnames(sampleDistMatrix) <- colnames(combat_edata)
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)
levels(sampleTable$species) = c('A.echinatior','C.biroi','D.quadriceps','L.humile','L.niger','M.pharaonis','S.invicta')
levels(sampleTable$caste) = c("Gyne",'Small worker','Non-reproductive','Reproductive','Worker')
sampleTable_col = data.frame(Species = sampleTable$species,row.names = row.names(sampleTable))
sampleTable_row = data.frame(Caste = sampleTable$caste,row.names = row.names(sampleTable))
sp_color = grey.colors(7,start = 0.1,end = 1)
caste_colour = rainbow(5)
ann_colors = list(Species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],M.pharaonis = sp_color[3],L.niger = sp_color[4],L.humile = sp_color[5],
                              C.biroi = sp_color[6],D.quadriceps = sp_color[7]),
                  Caste = c(Gyne = rgb(1,0,0,0.8),`Small worker` = 'yellow',Worker = rgb(0,0,1,0.8),
                            Reproductive = rgb(1,0,0,0.5), `Non-reproductive` = rgb(0,0,1,0.5)))

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=as.dist(1-cor(combat_edata,method = 's')),
         clustering_distance_cols=as.dist(1-cor(combat_edata,method = 's')),
         col=colors,annotation_colors = ann_colors)

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=dist(t(combat_edata)),
         clustering_distance_cols=dist(t(combat_edata)),
         col=colors,annotation_colors = ann_colors)

library(ggplot2)
se <- SummarizedExperiment(combat_edata - rowMeans(combat_edata),colData=sampleTable)
pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 7266, returnData=TRUE)

pcaData$species = factor(pcaData$species,levels = c('A.echinatior',"S.invicta","M.pharaonis","L.niger",'L.humile','C.biroi','D.quadriceps'))
pcaData$caste = factor(pcaData$caste,levels = c("Gyne","Small worker",'Worker', "Reproductive","Non-reproductive"))
levels(pcaData$caste) = c("Gyne","Worker","Worker","Reproductive","Non-reproductive")

percentVar_all <- round(100 * attr(pcaData, "percentVar"))
names(pcaData)[c(4,5)] = c('Species',"Caste")
ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
  geom_point(size=3,alpha = 1) +
  scale_shape_manual(values = c(15:17,0,1,2,5))+
  scale_color_manual(values = c('red','blue','purple','black'))+
  xlab(paste0("PC1 (",percentVar_all[1],"%)")) +
  ylab(paste0("PC2 (",percentVar_all[2],"%)")) + 
  coord_fixed()
