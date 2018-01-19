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
#####
ortholog_exp = cbind(Aech_exp[match(gene_ortholog_table$Aech, rownames(Aech_exp)),],
                     Lhum_exp[match(gene_ortholog_table$Lhum, rownames(Lhum_exp)),],
                     Mpha_exp[match(gene_ortholog_table$Mpha, rownames(Mpha_exp)),],
                     Sinv_exp[match(gene_ortholog_table$Sinv, rownames(Sinv_exp)),],
                     Lnig_exp[match(gene_ortholog_table$Lnig, rownames(Lnig_exp)),])

caste = factor(c(rep(c('gyne','minor'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4),rep(c('gyne','worker'),4)))

colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),
                  rep(c(16:19),each = 2),rep(c(20:23),each = 2)))


species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8),rep("Lnig",8)))


sampleTable <- data.frame(caste = caste, species = species_info, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(ortholog_exp)

library("pheatmap")
#ortholog_exp = ortholog_exp[apply(ortholog_exp, 1, function(x) length(x[x> .05])>=3),]
exp = log2(as.matrix(ortholog_exp) + 1e-5)
exp = normalize.quantiles(exp)

row.names(exp) = row.names(ortholog_exp)
colnames(exp) = colnames(ortholog_exp)
exp = exp[!apply(exp, 1, anyNA),]

library("RColorBrewer")

# Try to combat to remove the species effect, then the colony effect.
sampleTable_T = sampleTable
exp_T = exp
library(sva)
#ortholog_exp_filtered = ortholog_exp[!apply(ortholog_exp,1, anyNA),] 
#subset_test = c(1:8)
subset_test_T = grep('Aech',rownames(sampleTable_T))
subset_test_T = grep('Sinv',rownames(sampleTable_T))
subset_test_T = grep('Mpha',rownames(sampleTable_T))
subset_test_T = grep('Lnig',rownames(sampleTable_T))
subset_test_T = grep('Lhum',rownames(sampleTable_T))
subset_test_T = 99

filter_table_T = sampleTable_T[-subset_test_T,]
#batch = droplevels(filter_table_T$species)
batch = droplevels(filter_table_T$colony)

modcombat = model.matrix(~1, data=filter_table_T)
combat_edata_train = ComBat(dat=exp_T[,-subset_test_T], batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)
combat_edata_train = combat_edata_train[!apply(combat_edata_train,1, anyNA),] 
sampleDists = cor(combat_edata_train,method = 's')
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(combat_edata_train)
colnames(sampleDistMatrix) <- colnames(combat_edata_train)
colors <- colorRampPalette( brewer.pal(9, "Blues"))(255)
levels(sampleTable_T$species) = c('A.echinatior','L.humile','L.niger','M.pharaonis','S.invicta')
levels(sampleTable_T$caste) = c("Gyne",'Small worker','Worker')
sampleTable_col = data.frame(Species = filter_table_T$species,row.names = row.names(filter_table_T))
sampleTable_row = data.frame(Caste = filter_table_T$caste,row.names = row.names(filter_table_T))
sp_color = grey.colors(7,start = 0.1,end = 1)
caste_colour = rainbow(5)
ann_colors = list(Species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],M.pharaonis = sp_color[3],L.niger = sp_color[4],
                              L.humile = sp_color[5],C.biroi = sp_color[6],D.quadriceps = sp_color[7]),
                  Caste = c(Gyne = rgb(1,0,0,0.8),`Small worker` = 'yellow',Worker = rgb(0,0,1,0.8),
                            Reproductive = rgb(1,0,0,0.5), `Non-reproductive` = rgb(0,0,1,0.5)))

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=as.dist(1-cor(combat_edata_train,method = 's')),#sampleDists,
         clustering_distance_cols=as.dist(1-cor(combat_edata_train,method = 's')),#sampleDists,
         col=colors)#,annotation_colors = ann_colors)

pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_colnames = F,show_rownames = F,
         clustering_distance_rows=dist(t(combat_edata_train)),#sampleDists,
         clustering_distance_cols=dist(t(combat_edata_train)),#sampleDists,
         col=colors)#,annotation_colors = ann_colors)

library(ggplot2)
se <- SummarizedExperiment(combat_edata_train - rowMeans(combat_edata_train),colData=filter_table_T)
pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 7266, returnData=TRUE)

#pcaData$species = factor(pcaData$species,levels = c('A.echinatior',"S.invicta","M.pharaonis",'L.niger','L.humile'))
percentVar <- round(100 * attr(pcaData, "percentVar"))
names(pcaData)[c(4,5)] = c('Species',"Caste")
ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
  geom_point(size=3,alpha = .6) +
  scale_shape_manual(values = c(15:17,1,9,10))+
  scale_color_manual(values = c('red','brown','blue','purple','darkblue'))+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) + 
  coord_fixed()

#####

