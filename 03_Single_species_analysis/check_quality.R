library("RColorBrewer")
library("pheatmap")
library(readr)
gene_count_matrix <- read_csv("deg_ballgown/Aech_ERCC//gene_count_matrix.csv") #Output of Stringtie (gene expression matrix with ERCC)
count_matrix = gene_count_matrix[,c(paste('Aech',c(1:2,'3x',4:9,'10x',11,'12x','13x',14,15),sep ='_'))]
count_matrix = as.matrix(count_matrix)
rownames(count_matrix) = gene_count_matrix$X1
ercc_gene = grep('ERCC',rownames(count_matrix) )
#####
sampleDists = as.dist(1-cor(count_matrix[ercc_gene,],method = 's')) #Examine the Spearman correlation coefficient of ERCC among samples.
#sampleDists = as.dist(1-cor(count_matrix[-ercc_gene,],method = 's')) #Examine the Spearman correlation coefficient of all genes among samples.
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
