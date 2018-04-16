library(devtools)
library(Biobase)
library(preprocessCore)
library(tximport)
library('DESeq2')
library("RColorBrewer")
library("pheatmap")
gene_ortholog_table = read.table('../ortholog_seven/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua'))

read_input = function(species_header, col_name){
  aech_files = c(paste(species_header,col_name,sep ='_'))
  files <- file.path("~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/deg_salmon_gemoma/", species_header,
                     'quants', aech_files, "quant.sf")
  names(files) <- aech_files
  tx2gene <- read.csv(paste('~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/deg_salmon_gemoma/', species_header,'/',species_header,
                            '_gemoma_t2g.txt',sep = ''),header = F, sep = '\t')
  tx2gene$V1 = toupper(tx2gene$V1)
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                         countsFromAbundance = 'lengthScaledTPM')
  abundance_data = txi.salmon$abundance
  rownames(abundance_data)  = paste(species_header,  rownames(abundance_data), sep = '_')
  return(abundance_data)
}

Aech_exp = read_input('Aech',col_name = c(4,6,7,9,'10x','12x','13x',15))
Lhum_exp = read_input('Lhum',col_name = c(3:10))
Mpha_exp = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10))
Sinv_exp = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16))
Lnig_exp = read_input('Lnig',col_name = c(1:8))
Cbir_exp = read_input('Cbir',col_name = c(1:8))
Dqua_exp = read_input('Dqua',col_name = c(1:5,7:13))
#####
ortholog_exp = cbind(Aech_exp[match(gene_ortholog_table$Aech, rownames(Aech_exp)),],
                     Lhum_exp[match(gene_ortholog_table$Lhum, rownames(Lhum_exp)),],
                     Mpha_exp[match(gene_ortholog_table$Mpha, rownames(Mpha_exp)),],
                     Sinv_exp[match(gene_ortholog_table$Sinv, rownames(Sinv_exp)),],
                     Lnig_exp[match(gene_ortholog_table$Lnig, rownames(Lnig_exp)),],
                     Cbir_exp[match(gene_ortholog_table$Cbir, rownames(Cbir_exp)),],
                     Dqua_exp[match(gene_ortholog_table$Dqua, rownames(Dqua_exp)),])
                     
caste = factor(c(rep(c('gyne','minor'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4),rep(c('gyne','worker'),4),
                 rep(c('repro','non-repro'),4) ,rep(c('repro','non-repro'), each = 6)))

colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),rep(c(16:19),each = 2),rep(c(16:19)+20,each = 2),
                  rep(c(20:23),each = 2), rep(c(24:29),2)))


species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8), rep("Lnig",8),rep('Cbir',8),rep('Dqua',12)))


sampleTable <- data.frame(caste = caste, species = species_info, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(ortholog_exp)

library("pheatmap")
#ortholog_exp = ortholog_exp[apply(ortholog_exp, 1, function(x) length(x[x> .05])>=3),]
exp_data = log2(ortholog_exp + 1e-5)

exp = normalize.quantiles(as.matrix(exp_data))
row.names(exp) = row.names(exp_data)
colnames(exp) = colnames(exp_data)
exp = exp[!apply(exp, 1, anyNA),]


# Try to combat to remove the species effect, then the colony effect.

library(sva)
ortholog_exp_filtered = ortholog_exp[!apply(ortholog_exp,1, anyNA),] 
subset_test = c(1:62)
subset_test = c(1:42)
filter_table = sampleTable[subset_test,]
#filter_table$batch_t = as.factor(c(rep(1,26),rep(2,36)))
#batch = droplevels(filter_table$species)
batch = droplevels(filter_table$colony)

modcombat = model.matrix(~1, data=filter_table)
combat_edata = ComBat(dat=exp[,subset_test], batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)

sampleDists = cor(combat_edata,method = 's')
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(combat_edata)
colnames(sampleDistMatrix) <- colnames(combat_edata)
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)
levels(sampleTable$species) = c('A.echinatior','C.biroi','D.quadriceps','L.humile','L.niger','M.pharaonis','S.invicta')
levels(sampleTable$caste) = c("Gyne",'Small worker','Non-reproductive','Reproductive','Worker')
sampleTable_col = data.frame(Species = filter_table$species,row.names = row.names(filter_table))
sampleTable_row = data.frame(Caste = filter_table$caste,row.names = row.names(filter_table))
sp_color = grey.colors(7,start = 0.1,end = 1)
caste_colour = rainbow(5)
ann_colors = list(Species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],M.pharaonis = sp_color[3],L.niger = sp_color[4],L.humile = sp_color[5],
                              C.biroi = sp_color[6],D.quadriceps = sp_color[7]),
                  Caste = c(Gyne = rgb(1,0,0,0.8),`Small worker` = 'yellow',Worker = rgb(0,0,1,0.8),
                            Reproductive = rgb(1,0,0,0.5), `Non-reproductive` = rgb(0,0,1,0.5)))

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=as.dist(1-cor(combat_edata,method = 's')),#sampleDists,
         clustering_distance_cols=as.dist(1-cor(combat_edata,method = 's')),#sampleDists,
         col=colors)#,annotation_colors = ann_colors)

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=dist(t(combat_edata)),#sampleDists,
         clustering_distance_cols=dist(t(combat_edata)),#sampleDists,
         col=colors)#,annotation_colors = ann_colors)

library(ggplot2)
se <- SummarizedExperiment(combat_edata - rowMeans(combat_edata),colData=filter_table)
pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 7266, returnData=TRUE)

pcaData$species = factor(pcaData$species,levels = c('A.echinatior',"S.invicta","M.pharaonis","L.niger",'L.humile','C.biroi','D.quadriceps'))
pcaData$caste = factor(pcaData$caste,levels = c("Gyne","Small worker",'Worker', "Reproductive","Non-reproductive"))

#percentVar <- round(100 * attr(pcaData, "percentVar"))
names(pcaData)[c(4,5)] = c('Species',"Caste")
ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
  geom_point(size=3,alpha = .6) +
  scale_shape_manual(values = c(15:17,0,1,9,10))+
  scale_color_manual(values = c('red','brown','blue','purple','darkblue'))+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) + 
  coord_fixed()

#####

