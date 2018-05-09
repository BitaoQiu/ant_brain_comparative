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
combat_edata = ComBat(dat=exp, batch=batch, mod=modcombat,mean.only = F,
                            par.prior=TRUE,  prior.plots=FALSE)
combat_edata = combat_edata[!apply(combat_edata,1, anyNA),]
combat_edata_net = cor(t(combat_edata),method = 's') #Using the output of colony-level transcriptome normalization, constructed Spearman correlation matrix

grn = read.table('input/candidate_deg.txt',header = T,sep = '\t') #Imported genes significantly differentially expressed between castes across ants
res_all = read.table('input/res_five_ants.tsv',header = T, stringsAsFactors = F) # Imported the DEGs caste-biased level
gene_age = read.table('input/Aech_age.txt',header = T,stringsAsFactors = F) # Imported gene ages.
Aech_t2g=  read.csv('input/deg_salmon_gemoma/Aech/Aech_gemoma_t2g.txt',header = F, sep = '\t', stringsAsFactors = F)
Aech_t2g$V1 = toupper(Aech_t2g$V1)
Aech_t2g$V2 = paste("Aech",Aech_t2g$V2, sep = '_')
gene_age$geneID = Aech_t2g$V2[match(rownames(gene_age), Aech_t2g$V1)]
gene_age = gene_age[which(gene_age$geneID %in% rownames(combat_edata_net)),]
gene_age$age = factor(gene_age$age,
                      levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
levels(gene_age$age) = c('Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Species specific')
gene_age$age[gene_age$age %in% c('Myrmicinae','Species specific')] = 'Formicidae'
filter_age = gene_age$geneID[which(gene_age$age %in% c('Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Species specific'))] #Filtered genes with deputed age
combat_edata_net.f = combat_edata_net[filter_age,filter_age]
#####
net_background = apply(abs(combat_edata_net),1,median)
net_connectivity.test = apply(combat_edata_net, 1,FUN = function(x){
  wilcox.test(abs(x),net_background,alternative = 'g',paired = T)$p.value}) # One-sided wilcox test to compare the connectivity of target genes with that of whole genome background.
net_connectivity.test.adj = p.adjust(net_connectivity.test)
net_connectivity.test.adj[which(net_connectivity.test.adj == min(net_connectivity.test.adj))]
#####
n = 500 #Randomly sampled 500 genes
## May also add log2 expression ratio for all the genes.
tmp_gene = c(rownames(grn), sample(rownames(combat_edata_net.f),n))
tmp_gene = sample(rownames(combat_edata_net.f),n)
caste_bias = data.frame(caste_bias_2 = res_all$log2FoldChange[match(tmp_gene,rownames(res_all))],
                        conservation = -log10(net_connectivity.test.adj[tmp_gene]),
                        row.names = tmp_gene)
caste_bias$caste_bias_2[which(caste_bias$caste_bias_2 == min(caste_bias$caste_bias_2))] = -max(caste_bias$caste_bias_2)
caste_bias$conservation[is.infinite(caste_bias$conservation)] = 300 
age_bias = data.frame(Evolutionary.origin = gene_age$age[match(tmp_gene, gene_age$geneID)],
                      row.names = tmp_gene)

tmp_col = gray.colors(8, start = 0, end = 1,gamma = 2.8)
color.ramp <- colorRamp(c("blue", "red"))( (0:4)/4 )

ann_colors = list( Evolutionary.origin = c(Neoptera =tmp_col[1],Endopterygota = tmp_col[2],Hymenoptera=tmp_col[3],Aprocrita = tmp_col[4],Aculeata = tmp_col[5],Formicidae = tmp_col[6]),
                   caste_bias_2 = colorRampPalette(c('blue','white',"red"),bias =1, interpolate = "linear")( 5 ))
pheatmap(abs(combat_edata_net.f[tmp_gene,tmp_gene]), cluster_rows=T, cluster_cols=T, show_rownames=F,show_colnames = F, 
         annotation_col = caste_bias, annotation_row = age_bias,annotation_colors =  ann_colors, 
         clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation')
