library("RColorBrewer")
quick_check = function(pathway_edata,sampleTable, GO_term_tmp,gene_set){
  path_tmp = pathway_edata[rownames(pathway_edata) %in% gene_set[GO_term_tmp][[1]],]
  sampleDists = cor(path_tmp,method = 's')
  sampleTable_row <- data.frame(caste = sampleTable$caste,row.names =  colnames(path_tmp))
  sampleTable_col <- data.frame(species = sampleTable$species,row.names =  colnames(path_tmp))
  levels(sampleTable_col$species) = c('A.echinatior','L.humile',"M.pharaonis","S.invicta")
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
  sp_color = grey.colors(4,start = 0.1,end = 1)
  ann_colors = list(species = c(A.echinatior = sp_color[1],S.invicta = sp_color[2],L.humile = sp_color[3],M.pharaonis = sp_color[4]),
                    caste = c(gyne = rgb(1,0,0,0.5),worker = rgb(0,0,1,0.5)))
  pheatmap(sampleDistMatrix,annotation_col = sampleTable_col,annotation_row = sampleTable_row,show_rownames = F,show_colnames = F,
           clustering_distance_rows = as.dist(1-cor(path_tmp,method = 's')), annotation_colors = ann_colors,
           clustering_distance_cols = as.dist(1-cor(path_tmp,method = 's')),color = colors,main = paste(GO_term_tmp,", N = ",dim(path_tmp)[1],sep = ''))
}
quick_check = function(pathway_edata,sampleTable, GO_term_tmp,gene_set){
  path_tmp = pathway_edata[rownames(pathway_edata) %in% gene_set[GO_term_tmp][[1]],]
  se <- SummarizedExperiment(path_tmp - rowMeans(path_tmp),colData=sampleTable)
  pcaData <- plotPCA(DESeqTransform( se ), intgroup=c("species", "caste"),ntop = 7266, returnData=TRUE)
  levels(pcaData$species) = c('A.echinatior','L.humile',"M.pharaonis","S.invicta")
  levels(pcaData$caste) = c("Gyne",'Worker')
  pcaData$species = factor(pcaData$species,levels = c('A.echinatior',"S.invicta","M.pharaonis",'L.humile'))
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  names(pcaData)[c(4,5)] = c('Species',"Caste")
  ggplot(pcaData, aes(PC1, PC2, color=Caste, shape=Species)) +
    geom_point(size=3,alpha = .5) +
    scale_shape_manual(values = c(14:17))+
    scale_color_manual(values = c('red','blue'))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle(paste(GO_term_tmp,", N = ",dim(path_tmp)[1],sep = ''))+
    coord_fixed()
}

library(GO.db)

getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
  ans <- ans[!ans %in% goids]
}

level1_BP_terms <- getAllBPChildren("GO:0008150")     # 23 terms
level2_BP_terms <- getAllBPChildren(level1_BP_terms)  # 256 terms
level3_BP_terms <- getAllBPChildren(level2_BP_terms)  # 3059 terms
level4_BP_terms <- getAllBPChildren(level3_BP_terms)  # 9135 terms
level5_BP_terms <- getAllBPChildren(level4_BP_terms)  # 15023 terms
level6_BP_terms <- getAllBPChildren(level5_BP_terms)  # 15023 terms
level7_BP_terms <- getAllBPChildren(level6_BP_terms)  # 15023 terms

select_BP = go_p_value_filter[which(strtrim(row.names(go_p_value_filter),10) %in% level1_BP_terms),]
select_BP[order(select_BP$cor,decreasing = T)[1:30],]
path_tmp_table
quick_check(pathway_edata,sampleTable,"GO:0040007 growth",gene_set = go.bp)


####
myIds <- strtrim(row.names(go_p_value_filter[which(go_p_value_filter$cor > 0.8 & go_p_value_filter$pvalue < 0.01),]),10)
myCollection <- GOCollection(myIds)
fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
GO_slim_out = goSlim(myCollection, slim, "BP")
GO_slim_out$Term = as.character(GO_slim_out$Term)
GO_slim_out = GO_slim_out[order(GO_slim_out$Count,decreasing = T)[2:10],]
GO_slim_out$Term <-  factor(GO_slim_out$Term, levels=GO_slim_out[order(GO_slim_out$Count,decreasing=TRUE),'Term'])
ggplot(data =GO_slim_out, aes(Term,Count) )+
  geom_bar(stat = 'identity')

