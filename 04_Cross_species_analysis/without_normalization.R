library(devtools)
library(Biobase)
library(preprocessCore)
library(tximport)
library('DESeq2')
library("RColorBrewer")
library("pheatmap")

gene_ortholog_table = read.table('input/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua')) #Read in the ortholog table (from Step 1)

read_input = function(species_header, col_name){
  aech_files = c(paste(species_header,col_name,sep ='_'))
  files <- file.path("input/deg_salmon_gemoma/", species_header,
                     'quants', aech_files, "quant.sf") # Output of salmon, transcript level quantification
  names(files) <- aech_files
  tx2gene <- read.csv(paste('input/deg_salmon_gemoma/', species_header,'/',species_header,
                            '_gemoma_t2g.txt',sep = ''),header = F, sep = '\t') # Transcript to gene information, for gene level quantification
  tx2gene$V1 = toupper(tx2gene$V1)
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                         countsFromAbundance = 'lengthScaledTPM')
  abundance_data = txi.salmon$abundance
  rownames(abundance_data)  = paste(species_header,  rownames(abundance_data), sep = '_')
  return(abundance_data)
}

Aech_exp = read_input('Aech',col_name = c(4,6,7,9,'10x','12x','13x',15)) #Importing small worker and gyne brain samples of A.echinatior (we removed samples 1:3 because ERCC suggested poor quality of these samples)
Lhum_exp = read_input('Lhum',col_name = c(3:10)) #Samples of L.humile (We removed sample 1 and 2 because this colony was collected from different area)
Mpha_exp = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10)) #Samples of M. pharaonis
Sinv_exp = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16)) #Samples of S.invicta
Lnig_exp = read_input('Lnig',col_name = c(1:8)) #Samples of L.niger
Cbir_exp = read_input('Cbir',col_name = c(1:8)) # Samples of O.biroi
Dqua_exp = read_input('Dqua',col_name = c(1:5,7:13)) # D.quadriceps, removed 2 samples without colony information
#####
#Expression matrix for all 1-to-1 orthologous genes (only for 5 typical ant species)
ortholog_exp = cbind(Aech_exp[match(gene_ortholog_table$Aech, rownames(Aech_exp)),],
                     Lhum_exp[match(gene_ortholog_table$Lhum, rownames(Lhum_exp)),],
                     Mpha_exp[match(gene_ortholog_table$Mpha, rownames(Mpha_exp)),],
                     Sinv_exp[match(gene_ortholog_table$Sinv, rownames(Sinv_exp)),],
                     Lnig_exp[match(gene_ortholog_table$Lnig, rownames(Lnig_exp)),])
                     #,Cbir_exp[match(gene_ortholog_table$Cbir, rownames(Cbir_exp)),],
                     #Dqua_exp[match(gene_ortholog_table$Dqua, rownames(Dqua_exp)),])
#Phenotype information: Gyne, Minor worker, Worker, Reproductive worker, and Non-reproductive worker
caste = factor(c(rep(c('gyne','minor'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4),rep(c('gyne','worker'),4)))
               #,rep(c('repro','non-repro'),4) ,rep(c('repro','non-repro'), each = 6)))
#Colony information.
colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),rep(c(16:19),each = 2),rep(c(16:19)+20,each = 2)))
                  #,rep(c(20:23),each = 2), rep(c(24:29),2)))

#Species information.
species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8), rep("Lnig",8)))
                       #,rep('Cbir',8),rep('Dqua',12)))

#Construct sample information table.
sampleTable <- data.frame(caste = caste, species = species_info, colony = colony)
rownames(sampleTable) <- colnames(ortholog_exp)

exp_data = log2(ortholog_exp + 1e-5) # Log transformation, add 1e-5 to avoid log0
exp = normalize.quantiles(as.matrix(exp_data)) # Quantile normalization
row.names(exp) = row.names(exp_data)
colnames(exp) = colnames(exp_data)
exp = exp[!apply(exp, 1, anyNA),] #Removed genes showing NA (e.g. without expression)


sampleDists = cor(exp,method = 's')

sampleTable_row <- data.frame(Caste = caste,row.names =  colnames(ortholog_exp))
sampleTable_col <- data.frame(Species = species_info,row.names =  colnames(ortholog_exp))
levels(sampleTable_col$Species) = c('A. echinatior','L. humile','L. niger',"M. pharaonis","S. invicta")
levels(sampleTable_row$Caste) = c("Gyne","Worker",'Worker') #Combine minor worker and worker as one caste, for the sake of presentation

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
sp_color = grey.colors(5,start = 0,end = 1)

ann_colors = list(Species = c('A. echinatior' = sp_color[1],'S. invicta' = sp_color[2],'M. pharaonis' = sp_color[3],'L. niger' = sp_color[4],'L. humile' = sp_color[5]),
                  Caste = c(Gyne = rgb(1,0,0,0.6),Worker = rgb(0,0,1,0.6)))

#Reordered the hclust output, so that within branch, gyne will be plotted first then workers, and then species, in order to avoid visual confusion. Note that the results are topologically the same before and after the reordering.
callback = function(hc, mat){
  caste = factor(c(rep(c('g','w'),4),rep(c('g','w'),4),rep(c('g','w'),5),rep(c('g','w'),4),rep(c('g','w'),4)))
  Species =  factor(c(rep("A",8),rep("Lu",8),rep("Mp",10),rep("Si",8),rep("Ln",8)))
  combine = factor(c(paste(caste,Species,sep = '.')))
  combine = factor(combine, levels = paste(rep(c('g','w'),5),rep(c('A','Mp','Si','Lu','Ln'),each = 2),sep = '.'))
  dend = reorder(as.dendrogram(hc), wts = combine)
  as.hclust(dend)
}

pheatmap(sampleDists,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_rownames = F,show_colnames = F, clustering_callback = callback,
         clustering_distance_rows = dist(t(exp)), annotation_colors = ann_colors,
         clustering_distance_cols=dist(t(exp)),color = colors) #Clustering using euclidean distance

pheatmap(sampleDists,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_rownames = F,show_colnames = F,legend  = T, clustering_callback = callback,
         clustering_distance_rows = as.dist(1-cor(exp,method = 's')), annotation_colors = ann_colors,
         clustering_distance_cols=as.dist(1-cor(exp,method = 's')),color = colors) #Clustering using correlation coefficient

library(ggplot2)
se <- SummarizedExperiment(exp,colData=sampleTable)
pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 10000, returnData=TRUE)
levels(pcaData$species) = c('A. echinatior','L. humile','L. niger',"M. pharaonis","S. invicta")
pcaData$species = factor(pcaData$species,levels = c('A. echinatior',"S. invicta","M. pharaonis",'L. niger','L. humile'))
names(pcaData)[c(4,5)] = c('Species','Caste')
levels(pcaData$Caste) = c("Gyne",'Worker','Worker')
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Plotting the output of PCA
ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
  geom_point(size=3,alpha = .8) +
  scale_shape_manual(values = c(15:17,0,1))+
  scale_color_manual(values = c('red','blue'))+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  coord_fixed()

