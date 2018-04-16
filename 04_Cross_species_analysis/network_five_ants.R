combat_edata_net = cor(t(combat_edata),method = 's') #Using the output of colony-level transcriptome normalization, constructed Spearman correlation matrix

library("pheatmap")
library("RColorBrewer")
grn = read.table('candidate_deg.txt',header = T,sep = '\t') #Imported genes significantly differentially expressed between castes across ants
res_all = read.table('res_five_ants.tsv',header = T, stringsAsFactors = F) # Imported the DEGs caste-biased level
gene_age = read.table('Aech_age.txt',header = T,stringsAsFactors = F) # Imported gene ages.
Aech_t2g=  read.csv('deg_salmon_gemoma/Aech/Aech_gemoma_t2g.txt',header = F, sep = '\t', stringsAsFactors = F)
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
