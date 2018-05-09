library(sva)
library(devtools)
library(Biobase)
library(preprocessCore)
library(tximport)
library('DESeq2')
library("RColorBrewer")
library("pheatmap")

gene_ortholog_table = read.table('input/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua')) #orthologous genes relations among seven ant species

read_input = function(species_header, col_name){
  aech_files = c(paste(species_header,col_name,sep ='_'))
  files <- file.path("input/deg_salmon_gemoma/", species_header,#Output of salmon
                     'quants', aech_files, "quant.sf")
  names(files) <- aech_files
  tx2gene <- read.csv(paste('input/deg_salmon_gemoma/', species_header,'/',species_header,
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
# Gene expression matrix for 1-to-1 orthologous genes in five typical ant species
ortholog_exp = cbind(Aech_exp[match(gene_ortholog_table$Aech, rownames(Aech_exp)),],
                     Lhum_exp[match(gene_ortholog_table$Lhum, rownames(Lhum_exp)),],
                     Mpha_exp[match(gene_ortholog_table$Mpha, rownames(Mpha_exp)),],
                     Sinv_exp[match(gene_ortholog_table$Sinv, rownames(Sinv_exp)),],
                     Lnig_exp[match(gene_ortholog_table$Lnig, rownames(Lnig_exp)),])

caste = factor(c(rep(c('gyne','minor'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4),rep(c('gyne','worker'),4)))

colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),
                  rep(c(16:19),each = 2),rep(c(20:23),each = 2)))


species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8),rep("Lnig",8)))


sampleTable <- data.frame(caste = caste, species = species_info, colony = colony)
rownames(sampleTable) <- colnames(ortholog_exp)

library("pheatmap")
exp = log2(as.matrix(ortholog_exp) + 1e-5)
exp = normalize.quantiles(exp)

row.names(exp) = row.names(ortholog_exp)
colnames(exp) = colnames(ortholog_exp)
exp = exp[!apply(exp, 1, anyNA),]



# Try to combat to remove the species effect, then the colony effect.

ortholog_exp_filtered = ortholog_exp[!apply(ortholog_exp,1, anyNA),] 

filter_table = sampleTable
#batch = droplevels(filter_table$species) #Species level normalization 
batch = droplevels(filter_table$colony)  #Colony level normalization, note that normalization for colony always normalized for species effect.

modcombat = model.matrix(~1, data=filter_table)
combat_edata_train = ComBat(dat=exp, batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)
combat_edata_train = combat_edata_train[!apply(combat_edata_train,1, anyNA),] #This is the normalized transcriptome from five ant species. 
