# Continue from the output of salmon_reference.R
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
Cbir_exp = read_input('Cbir',col_name = c(1:8))
Dqua_exp = read_input('Dqua',col_name = c(1:5,7:13))
#####
# Gene expression matrix for 1-to-1 orthologous genes in all seven ant species
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
exp_data = log2(ortholog_exp + 1e-5)
exp = normalize.quantiles(as.matrix(exp_data))
row.names(exp) = row.names(exp_data)
colnames(exp) = colnames(exp_data)
exp = exp[!apply(exp, 1, anyNA),]


library(sva)
ortholog_exp_filtered = ortholog_exp[!apply(ortholog_exp,1, anyNA),] 

filter_table = sampleTable
#batch = droplevels(filter_table$species)
batch = droplevels(filter_table$colony) #Note that should choose the same level of normalization for both reference and target.

modcombat = model.matrix(~1, data=filter_table)
combat_edata_target = ComBat(dat=exp, batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE) #Normalized expression level of both reference and target species.
save(combat_edata_target, filter_table, file = "03_input_reference_2.RData")
