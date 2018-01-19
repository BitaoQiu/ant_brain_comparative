#gene_age
t2gene = read.table('~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/ortholog/ortholog_gexpr_TPM_four_sp/candidate_gene/Aech_gemoma_t2g.txt')
library(readr)
tblast <- read_delim("~/Library/Mobile Documents/com~apple~CloudDocs/KU_data_analysis/gene_age/gene_age_v2/ortholog/ortholog_gexpr_TPM_four_sp/candidate_gene/proteins.fasta.blastout", 
                     "\t", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
t2gene$V1 = toupper(t2gene$V1)
t2gene$V2 = paste("Aech",t2gene$V2,sep = '_')
head(t2gene)
head(tblast)
t2gene$blast = tblast$X17[match(t2gene$V1,tblast$X1)]
t2gene$evalue = tblast$X15[match(t2gene$V1,tblast$X1)]
t2gene$coverage = tblast$X11[match(t2gene$V1,tblast$X1)]
library("ggplot2")
library(devtools)
library(Biobase)
library(preprocessCore)
library(tximport)
library('DESeq2')
library("RColorBrewer")
library("pheatmap")
gene_ortholog_table = read.table('../../ortholog/ortholog_seven/gene_table.poff', col.names = c('Aech','Sinv','Mpha','Lnig','Lhum','Cbir','Dqua'))

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
#Aech_exp = read_input('Aech',col_name = c(4,5,7,8,'10x',11,'13x',14))
Lhum_exp = read_input('Lhum',col_name = c(3:10))
#Lhum_exp = read_input('Lhum',col_name = c(1:10))
Mpha_exp = read_input('Mpha',col_name = c(1:6,'7x',8,'9x',10))
Sinv_exp = read_input('Sinv',col_name = c(3,4,7,8,11,12,15,16))
Lnig_exp = read_input('Lnig',col_name = c(1:8))
#####
Aech_m = match(gene_ortholog_table$Aech, rownames(Aech_exp$abundance))
Lhum_m = match(gene_ortholog_table$Lhum, rownames(Lhum_exp$abundance))
Mpha_m = match(gene_ortholog_table$Mpha, rownames(Mpha_exp$abundance))
Sinv_m = match(gene_ortholog_table$Sinv, rownames(Sinv_exp$abundance))
Lnig_m = match(gene_ortholog_table$Lnig, rownames(Lnig_exp$abundance))
ortholog_exp = list(abundance = cbind(Aech_exp$abundance[Aech_m,],
                                      Lhum_exp$abundance[Lhum_m,],
                                      Mpha_exp$abundance[Mpha_m,],
                                      Sinv_exp$abundance[Sinv_m,],
                                      Lnig_exp$abundance[Lnig_m,]),
                    counts = cbind(Aech_exp$counts[Aech_m,],
                                   Lhum_exp$counts[Lhum_m,],
                                   Mpha_exp$counts[Mpha_m,],
                                   Sinv_exp$counts[Sinv_m,],
                                   Lnig_exp$counts[Lnig_m,]),
                    length = cbind(Aech_exp$length[Aech_m,],
                                   Lhum_exp$length[Lhum_m,],
                                   Mpha_exp$length[Mpha_m,],
                                   Sinv_exp$length[Sinv_m,],
                                   Lnig_exp$length[Lnig_m,]),
                    countsFromAbundance = 'lengthScaledTPM')
ortholog_exp$length = ortholog_exp$length[!apply(ortholog_exp$abundance,1,anyNA),]
ortholog_exp$counts = ortholog_exp$counts[!apply(ortholog_exp$abundance,1,anyNA),]
ortholog_exp$abundance = ortholog_exp$abundance[!apply(ortholog_exp$abundance,1,anyNA),]

batch = factor(c(rep(c('gyne','worker'),4),rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4),rep(c('gyne','worker'),4)))
#batch = factor(c(rep(c('gyne','worker'),4),rep(c('gyne','worker'),5),rep(c('gyne','worker'),5),rep(c('gyne','worker'),4)))

colony = factor(c(rep(c(2:5),each = 2),rep(c(7:10),each = 2),rep(c(11:15),each = 2),
                  rep(c(16:19),each = 2),rep(c(20:23),each = 2)))
#colony.n = factor(c(rep(c(2:5),each = 2),rep(c(2:5),each = 2),rep(c(2:6),each = 2),
 #                   rep(c(2:5),each = 2)))
#colony = factor(c(rep(c(2:5),each = 2),rep(c(6:10),each = 2),rep(c(11:15),each = 2),
#                  rep(c(16:19),each = 2)))
species_info =  factor(c(rep("Aech",8),rep("Lhum",8),rep("Mpha",10),rep("Sinv",8),rep("Lnig",8)))
#species_info =  factor(c(rep("Aech",8),rep("Lhum",10),rep("Mpha",10),rep("Sinv",8)))

sampleTable <- data.frame(caste = batch, species = species_info, colony = colony)#,     type = c(rep('A',15),rep('B',4)))
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

Lnig_sampleTable <- sampleTable[c(35:42),]
Lnig_dds = DESeqDataSetFromTximport(Lnig_exp,Lnig_sampleTable, ~caste + colony)
Lnig_dds = estimateSizeFactors(Lnig_dds)

colData(Aech_dds)
#####
dds <- DESeqDataSetFromTximport(ortholog_exp, sampleTable, 
                                ~ colony + caste)# + type )

colData(dds)$sizeFactor = c(colData(Aech_dds)$sizeFactor,colData(Lhum_dds)$sizeFactor,colData(Mpha_dds)$sizeFactor,colData(Sinv_dds)$sizeFactor,colData(Lnig_dds)$sizeFactor)
dds = estimateDispersions(dds)
dds =  nbinomWaldTest(dds)

res <- results(dds,contrast = c('caste','gyne','worker'),alpha = 1e-2,lfcThreshold = log2(1.5))
summary(res)
#res_filter = res[which(res$padj < 0.01 | rownames(res) %in% c('Aech_gene_13482','Aech_gene_9986')),]
res_filter = res[which(res$padj <1e-2),]
#res_filter = res_filter[which(abs(res_filter$log2FoldChange)>1),]
res_filter$name = t2gene$blast[match(rownames(res_filter),t2gene$V2)]
res_filter = data.frame(res_filter)
res_filter
res_ver1 = data.frame(res_filter[order(res_filter$log2FoldChange),])
res_ver2 = data.frame(res_filter[order(res_filter$log2FoldChange),])
test_gene = c('Aech_gene_5241')
d <- plotCounts(dds,normalized = T, gene=test_gene, intgroup=c("caste",'species','colony'),returnData = T)
ggplot(data = d, aes(x = caste,y = count, colour = species,group = colony))+
  geom_point()+
  geom_line()+
  scale_y_log10()

make_gene_table = function(gene_list){
  out_d = data.frame(matrix(ncol = 6, nrow = 0))
  names(out_d) = c('count','caste','species','colony','geneID','annotation')
  for(gene in gene_list){
    tmp_d = plotCounts(dds,normalized = T, gene=gene, intgroup=c("caste",'species','colony'), 
                       returnData=TRUE)
    tmp_d$geneID = gene
    tmp_d$annotation = gsub('PREDICTED: ','',t2gene$blast[which(t2gene$V2 == gene)[1]])
    tmp_d$annotation = gsub("\\[Trachymyrmex septentrionalis]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Acromyrmex echinatior]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Wasmannia auropunctata]",'',tmp_d$annotation)
    tmp_d$annotation = gsub("\\[Trachymyrmex zeteki]",'',tmp_d$annotation)
    
    out_d = rbind(out_d,tmp_d)}
  return(out_d)
}
test_gene_list = c('Aech_gene_13364','Aech_gene_4170','Aech_gene_13934','Aech_gene_703',
                   'Aech_gene_15227','Aech_gene_5241','Aech_gene_9986','Aech_gene_13482',
                   'Aech_gene_5533','Aech_gene_3189','Aech_gene_5241')#,'Aech_gene_1056')
test_gene_list = c('Aech_gene_13364','Aech_gene_4170','Aech_gene_13934','Aech_gene_703',
                   'Aech_gene_15227','Aech_gene_5241','Aech_gene_9986','Aech_gene_13482','Aech_gene_5533','Aech_gene_3189','Aech_gene_840')#,'Aech_gene_1056')
res_filter[,c(2,6,7)]
write.table(res_filter[,c(2,6,7)],'../../ortholog/ortholog_gexpr_TPM_five_sp/candidate_deg.txt',quote = F,row.names = T, col.names = T, sep = '\t')
test_gene_list = c('Aech_gene_3189','Aech_gene_5319','Aech_gene_8046','Aech_gene_13364','Aech_gene_15227','Aech_gene_5533','Aech_gene_840','Aech_gene_5241')

d_test_gene_list = make_gene_table(test_gene_list)
names(d_test_gene_list)[c(2,3)] = c("Caste",'Species')
levels(d_test_gene_list$Species) = c("A.echinatior",'L.humile','L.niger','M.pharaonis','S.invicta')
d_test_gene_list$Species = factor(d_test_gene_list$Species, levels = c("A.echinatior",'M.pharaonis','S.invicta','L.niger','L.humile'))
levels(d_test_gene_list$Caste) = c('Gyne','Worker')
d_test_gene_list$fold = res$log2FoldChange[match(d_test_gene_list$geneID,rownames(res))]
d_test_gene_list$annotation = factor(d_test_gene_list$annotation)
d_test_gene_list$annotation = reorder(d_test_gene_list$annotation,X = -d_test_gene_list$fold)
levels(d_test_gene_list$annotation)[6] = 'IGFALS isoform X1'
#levels(d_test_gene_list$annotation)[7] = 'GABA receptor\n subunit beta isoform X3'

levels(d_test_gene_list$annotation)[7] = 'RERG-like protein'

ggplot(d_test_gene_list, aes(x=Caste, y=count,col = Species)) + 
  geom_jitter(height = 0,width = 0) + 
  geom_smooth(aes(group = Species),size=0.5,method = 'glm',se = F,linetype=1)+
  scale_y_log10()+
  theme(legend.position="top")+
  facet_wrap(~ annotation,scales = 'free_y',nrow = 2, dir="h",shrink = F)+
  ylab('Normalized count reads')

#Mpha_exp$abundance[rownames(Mpha_exp$abundance) == 'Mpha_gene_10588',]
#write.table(a,'~/Downloads/bias_four.txt',sep = '\t',quote = F,row.names = T)


rld <- rlog(dds, blind=FALSE)

vsd <- varianceStabilizingTransformation(dds, blind=T)
vsd.fast <- vst(dds, blind=T)
tmp = rld
tmp = vsd.fast
filter_table = sampleTable
batch = droplevels(filter_table$species)
batch = droplevels(filter_table$colony)
library(sva)
modcombat = model.matrix(~1, data=filter_table)
combat_edata = ComBat(dat=assay(tmp), batch=batch, mod=modcombat,mean.only = F,
                      par.prior=TRUE,  prior.plots=FALSE)


se <- SummarizedExperiment(combat_edata - rowMeans(combat_edata),colData=sampleTable)
pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 7266, returnData=TRUE)

levels(pcaData$species) = c("A.echinatior",'L.humile','L.niger','M.pharaonis','S.invicta')
pcaData$species = factor(pcaData$species,levels = c('A.echinatior','S.invicta',"M.pharaonis",'L.niger','L.humile'))
levels(pcaData$caste)  = c("Gyne",'Worker')
percentVar <- round(100 * attr(pcaData, "percentVar"))
names(pcaData)[c(4,5)] = c('Species',"Caste")
ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
  geom_point(size=3,alpha = .8) +
  scale_shape_manual(values = c(15:17,0,1))+
  scale_color_manual(values = c('red','blue','brown','blue','purple','darkblue'))+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) 
  