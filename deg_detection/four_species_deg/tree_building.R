library(phytools)
library(ape)
library(phangorn)
library(seqinr)
ant_tree = read.newick(file = 'ant_tree.nex')
names(stat_deg) = c("Aech","Aech",'Mpha','Sinv','Lhum','Cbir','Dqua')

plot(ant_tree)

dm = as.dist(1-cor(log2_tmp, use = 'c',method = 's'))
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)
layout(matrix(c(1,1), 1, 1), height=c(1,2))
par(mar = c(5,5,5,5)+ 0.1)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")
tmp_data = log2_deg[!apply(log2_deg ,1 , anyNA),]
tmp_data = tmp_data[,c(7,6,5,3,4,2)]
colnames(tmp_data) = c('Dqua','Cbir','Lhum','Mpha','Sinv','Aech')
tmp_dist = as.matrix(cor((tmp_data),use = 'c',method = 's'))
test = as.dist(1-cor((tmp_data),use = 'c',method = 's'))
plot(hclust(test))


test_data = t(log2_tmp)
test_data = test_data[,!apply(test_data,2, anyNA)]
phylosig(ant_tree, test_data[,1], method="lambda", test=T, se=NULL, start=NULL,
         control=list())


test_data_lambda  = apply(test_data,2,FUN = function(x){phylosig(ant_tree, x, method="lambda", test=T, nsim=1000, se=NULL, start=NULL,
                                            control=list())$lambda})
test_data_P  = apply(test_data,2,FUN = function(x){phylosig(ant_tree, x, method="lambda", test=T, nsim=1000, se=NULL, start=NULL,
                                                            control=list())$P})
PatristicDistMatrix<-cophenetic(ant_tree)
plot(PatristicDistMatrix)

tmp_data = log2_tmp[!apply(log2_tmp ,1 , anyNA),]
tmp_data = tmp_data[,c(6,5,4,2,3,1)]
colnames(tmp_data) = c('Dqua','Cbir','Lhum','Mpha','Sinv','Aech')
tmp_dist = as.matrix(cor((tmp_data),use = 'c',method = 's'))
plot(PatristicDistMatrix,tmp_dist, xlim = c(50,280),pch = 20, xlab = 'Phylogenetic distance(MYA)', ylab = 'spearman correlation of caste difference')
