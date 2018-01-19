library(devtools)
library(Biobase)
library(preprocessCore)
library(tximport)
library('DESeq2')
library("RColorBrewer")
library("pheatmap")
gene_ortholog_table = read.table('../../ortholog/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lhum','Cbir','Dqua'))

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
  #count_data = txi.salmon$counts
  rownames(txi.salmon$abundance)  = paste(species_header,  rownames(txi.salmon$abundance), sep = '_')
  rownames(txi.salmon$counts)  = paste(species_header,  rownames(txi.salmon$counts), sep = '_')
  rownames(txi.salmon$length)  = paste(species_header,  rownames(txi.salmon$length), sep = '_')
  return(txi.salmon)
}
Aech_exp = read_input('Aech',col_name = c(4,6,7,9,'10x','12x','13x',15))

Aech_exp = read_input('Aech',col_name = c(4,5,7,8,'10x',11,'13x',14))
Lhum_exp = read_input('Lhum',col_name = c(3:10))
#Lhum_exp = read_input('Lhum',col_name = c(1:10))
Mpha_exp = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10))
Sinv_exp = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16))
#####
Aech_m = match(gene_ortholog_table$Aech, rownames(Aech_exp$abundance))
Lhum_m = match(gene_ortholog_table$Lhum, rownames(Lhum_exp$abundance))
Mpha_m = match(gene_ortholog_table$Mpha, rownames(Mpha_exp$abundance))
Sinv_m = match(gene_ortholog_table$Sinv, rownames(Sinv_exp$abundance))
ortholog_exp = list(abundance = cbind(Aech_exp$abundance[Aech_m,],
                                      Lhum_exp$abundance[Lhum_m,],
                                      Mpha_exp$abundance[Mpha_m,],
                                      Sinv_exp$abundance[Sinv_m,]),
                     counts = cbind(Aech_exp$counts[Aech_m,],
                                   Lhum_exp$counts[Lhum_m,],
                                   Mpha_exp$counts[Mpha_m,],
                                   Sinv_exp$counts[Sinv_m,]),
                    length = cbind(Aech_exp$length[Aech_m,],
                                   Lhum_exp$length[Lhum_m,],
                                   Mpha_exp$length[Mpha_m,],
                                   Sinv_exp$length[Sinv_m,]),
                    countsFromAbundance = 'lengthScaledTPM')
ortholog_exp$length = ortholog_exp$length[!apply(ortholog_exp$abundance,1,anyNA),]
ortholog_exp$counts = ortholog_exp$counts[!apply(ortholog_exp$abundance,1,anyNA),]
ortholog_exp$abundance = ortholog_exp$abundance[!apply(ortholog_exp$abundance,1,anyNA),]

batch = factor(c(rep(c('gyne','worker'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4)))
#batch = factor(c(rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4)))

colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),
                  rep(c(16:19),each = 2)))
colony.n = factor(c(rep(c(2:5),each = 2),rep(c(2:5),each = 2),rep(c(2:6),each = 2),
                  rep(c(2:5),each = 2)))
#colony = factor(c(rep(c(2:5),each = 2),rep(c(6:10),each = 2),rep(c(11:15),each = 2),
#                  rep(c(16:19),each = 2)))
species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8)))
#species_info =  factor(c(rep("Aech",8),rep("Lhum",10),rep("Mpha",10),rep("Sinv",8)))

sampleTable <- data.frame(caste = batch, species = species_info, colony = colony,colony.n= colony.n)#,     type = c(rep('A',15),rep('B',4)))
rownames(sampleTable) <- colnames(ortholog_exp$abundance)

#####
Aech_sampleTable <- sampleTable[c(1:8),]
Aech_dds = DESeqDataSetFromTximport(Aech_exp,Aech_sampleTable, ~caste + colony)
Aech_dds = estimateSizeFactors(Aech_dds)

Lhum_sampleTable <- sampleTable[c(9:16),]
Lhum_dds = DESeqDataSetFromTximport(Lhum_exp,Lhum_sampleTable, ~caste + colony)
Lhum_dds = estimateSizeFactors(Lhum_dds)

Mpha_sampleTable <- sampleTable[c(17:26),]
Mpha_dds = DESeqDataSetFromTximport(Mpha_exp,Mpha_sampleTable, ~caste + colony)
Mpha_dds = estimateSizeFactors(Mpha_dds)

Sinv_sampleTable <- sampleTable[c(27:34),]
Sinv_dds = DESeqDataSetFromTximport(Sinv_exp,Sinv_sampleTable, ~caste + colony)
Sinv_dds = estimateSizeFactors(Sinv_dds)

colData(Aech_dds)
#####
dds <- DESeqDataSetFromTximport(ortholog_exp, sampleTable, 
                                ~caste + colony)# + type )

#dds <- DESeq(dds)
colData(dds)$sizeFactor = c(colData(Aech_dds)$sizeFactor,colData(Lhum_dds)$sizeFactor,colData(Mpha_dds)$sizeFactor,colData(Sinv_dds)$sizeFactor)
dds = estimateDispersions(dds)
dds =  nbinomWaldTest(dds)

res <- results(dds, c('caste','gyne','worker'),alpha = 0.005,lfcThreshold = log2(1.5))
summary(res)
res_filter = res[which(res$padj < 0.005 | rownames(res) %in% c('Aech_gene_13482','Aech_gene_9986')),]
#res_filter = res_filter[which(abs(res_filter$log2FoldChange)>1),]
res_filter$name = t2gene$blast[match(rownames(res_filter),t2gene$V2)]
res_filter = data.frame(res_filter)
res_ver3 = data.frame(res_filter[order(res_filter$log2FoldChange,decreasing = T),])
write.table(data.frame(res_filter),'../Amel/ant_res.txt',sep = '\t',quote = F,row.names = T, col.names = T)
res_LFC = lfcShrink(dds, contrast=c("caste","gyne","worker"), res=res)
summary(res)
res_a <- results(dds, contrast =list("caste_worker_vs_gyne") ,alpha = 0.01)
res_e = results(dds, contrast =list("caste_worker_vs_gyne","casteworker.speciesLhum") ,alpha = 0.01)
res_b <- results(dds, contrast =list("casteworker.speciesLhum") ,alpha = 0.01)
res_c <- results(dds, contrast =list("casteworker.speciesMpha") ,alpha = 0.01)
res_d <- results(dds, contrast =list("casteworker.speciesSinv") ,alpha = 0.01)
worker_bias = rbind(res_a[which(res_a$padj < 0.1 & res_a$log2FoldChange > 0 & res_b$padj > 0.2  &  res_c$padj > 0.2   & res_d$padj > 0.2),],
                    res_a[which(res_a$padj < 0.1 & res_a$log2FoldChange > 0 & res_b$log2FoldChange > 0 &  res_c$log2FoldChange > 0  & res_d$log2FoldChange > 0),])
worker_bias$name = t2gene$blast[match(rownames(worker_bias),t2gene$V2)]
gyen_bias = rbind(res_a[which(res_a$padj < 0.1 & res_a$log2FoldChange < 0 & res_b$padj > 0.2  &  res_c$padj > 0.2   & res_d$padj > 0.2),],
                  res_a[which(res_a$padj < 0.1 & res_a$log2FoldChange < 0 & res_b$log2FoldChange < 0 &  res_c$log2FoldChange < 0  & res_d$log2FoldChange < 0),])
gyen_bias$name = t2gene$blast[match(rownames(gyen_bias),t2gene$V2)]
data.frame(gyen_bias)
data.frame(worker_bias)


#plotMA(res, ylim=c(-2,2))
res_sig  = res[which(res$padj < 0.01 & abs(res$log2FoldChange)>1),]
res_sig  = res_LFC[which(res_LFC$padj < 0.01 & abs(res_LFC$log2FoldChange)>1),]
summary(res)
summary(res_sig)

table(abs(res_sig$log2FoldChange) >1)
res_sig$dmel = as.character(GO_annotation$name[match(rownames(res_sig),GO_annotation$V2)])
res_sig$id = t2gene$blast[match(rownames(res_sig),t2gene$V2)]
data.frame(res_sig[res_sig$log2FoldChange > 0,c(2,6,8)])
data.frame(res_sig[,c(2,5,7)])
a = data.frame(res_sig[,c(2,5,7)])
#res_sig  = res[which(res$padj < 0.0001 & abs(res$log2FoldChange)> 1),]
res_sig[order(res_sig$padj),]
res_LFC = lfcShrink(dds, contrast=c("caste","gyne","worker"), res=res)
plotMA(res_LFC, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup=c("caste"))

test_gene = c('Aech_gene_8046')
d <- plotCounts(dds,normalized = T, gene=test_gene, intgroup=c("caste",'species','colony'), 
                returnData=TRUE)
make_gene_table = function(gene_list){
  out_d = data.frame(matrix(ncol = 6, nrow = 0))
  names(out_d) = c('count','caste','species','colony','geneID','annotation')
  for(gene in gene_list){
    tmp_d = plotCounts(dds,normalized = T, gene=gene, intgroup=c("caste",'species','colony'), 
                            returnData=TRUE)
    tmp_d$geneID = gene
    tmp_d$annotation = gsub('PREDICTED: ','',t2gene$blast[which(t2gene$V2 == gene)][1])
    tmp_d$annotation = gsub("\\[Trachymyrmex septentrionalis]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Acromyrmex echinatior]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Wasmannia auropunctata]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Trachymyrmex zeteki]",'',tmp_d$annotation)
    
    out_d = rbind(out_d,tmp_d)}
  return(out_d)
}
test_gene_list = c('Aech_gene_13364','Aech_gene_4170','Aech_gene_13934','Aech_gene_703',
                   'Aech_gene_15227','Aech_gene_5241','Aech_gene_9986','Aech_gene_13482','Aech_gene_5533','Aech_gene_3189','Aech_gene_840')#,'Aech_gene_1056')
d_test_gene_list = make_gene_table(test_gene_list)
names(d_test_gene_list)[c(2,3)] = c("Caste",'Species')
levels(d_test_gene_list$Species) = c("A.echinatior",'L.humile','M.pharaonis','S.invicta')
d_test_gene_list$Species = factor(d_test_gene_list$Species, levels = c("A.echinatior",'M.pharaonis','S.invicta','L.humile'))
levels(d_test_gene_list$Caste) = c('Gyne','Worker')
d_test_gene_list$fold = res$log2FoldChange[match(d_test_gene_list$geneID,rownames(res))]
d_test_gene_list$annotation = factor(d_test_gene_list$annotation)
d_test_gene_list$annotation = reorder(d_test_gene_list$annotation,X = -d_test_gene_list$fold)
levels(d_test_gene_list$annotation)[2] = 'general OBP 69a-like'
levels(d_test_gene_list$annotation)[5] = 'GABA receptor subunit beta isoform X3'
levels(d_test_gene_list$annotation)[6] = 'CaM-kinase II'
levels(d_test_gene_list$annotation)[7] = 'RERG-like protein'
#d_test_gene_list$Bias = rep(c('Gyne-biased',"Worker-biased"),each = 136)
ggplot(d_test_gene_list, aes(x=Caste, y=count,col = Species)) + 
  geom_jitter(height = 0,width = 0) + 
  geom_smooth(aes(group = Species),size=1,method = 'glm',se = F,linetype=1)+
  scale_y_log10()+
  facet_wrap(~ annotation,scales = 'free_y',nrow = 3, dir="h",shrink = F)+
  ylab('Normalized count reads')

#Mpha_exp$abundance[rownames(Mpha_exp$abundance) == 'Mpha_gene_10588',]
#write.table(a,'~/Downloads/bias_four.txt',sep = '\t',quote = F,row.names = T)


rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("caste","species")])
pheatmap(assay(vsd.fast)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
plotPCA(vsd.fast, intgroup=c("caste", "species"))
