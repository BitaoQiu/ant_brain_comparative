w = read.table('dnds.tsv')
w = read.table('M0.txt')
w = read.table('M0_10.txt')
log2_deg$var = apply(log2_deg[,c(2:5)],1,sd)
log2_deg$w = w$V2[match(rownames(log2_deg),w$V1)]
log2_deg$fold = res_LFC$log2FoldChange[match(rownames(log2_deg),rownames(res_LFC))]
log2_deg$Exp = res_LFC$baseMean[match(rownames(log2_deg),rownames(res_LFC))]

log2_deg$caste = 'NDE'
log2_deg$caste[which(log2_deg$fold > 1)] = 'Gyne'
log2_deg$caste[which(log2_deg$fold < -1)] = 'Worker'
plot(log2_deg$var, log2_deg$w, col = c('red','black','blue')[(log2_deg$caste)],pch = 20,log = 'xy')
points(log2_deg$var, log2_deg$w, col = c('red',NA,'blue')[(log2_deg$caste)],pch = 20)
cor(log2_deg$var, log2_deg$fold,use = 'c',method = 's')
cor(log2_deg$w, log2_deg$fold,use = 'c',method = 's')
cor(log2_deg$var, log2_deg$w,use = 'c',method = 's')
cor(log2_deg$fold, log2_deg$Exp,use = 'c',method = 's')

#####
library(ggplot2)
library(MASS)
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
log2_deg_plot = log2_deg[!apply(log2_deg,1,anyNA),]
log2_deg_plot$density  = get_density(log2_deg_plot$var,log2_deg_plot$w)
log2_deg_plot = log2_deg_plot[-which(abs(log2_deg_plot$fold) < 1e-2),]
cor(log2_deg_plot$fold, log2_deg_plot$w,use = 'c',method = 's')

ggplot(data= log2_deg_plot,aes(x = fold, y = w, colour = caste))+
  geom_point()+
#  scale_color_distiller(palette = "RdBu")+
  xlab('Variation of caste-biased among ants')+
  geom_smooth(colour ='red',method = 'glm')+
  ylab('dN/dS')+
#  scale_x_sqrt()+
  scale_y_sqrt()

library(ppcor)
pcor.test(log2_deg_plot$var,log2_deg_plot$w,log(log2_deg_plot$Exp),method = 'p')

res$w = w$V2[match(rownames(res),w$V1)]
head(res)
w[w$V1 == 'Aech_gene_7060',]
deg_overlap_w$w = w$V2[match(rownames(deg_overlap_w),w$V1)]
deg_overlap_g$w = w$V2[match(rownames(deg_overlap_g),w$V1)]
deg_overlap_nde$w = w$V2[match(rownames(deg_overlap_nde),w$V1)]

plot(abs(res$log2FoldChange),res$w)
boxplot(deg_overlap_w$w,deg_overlap_g$w , deg_overlap_nde$w,names=c('Worker','Gyne','all'), main = 'One ratio model (three out of four degs)')
wilcox.test(deg_overlap_g$w, deg_overlap_nde$w)
wilcox.test(deg_overlap_w$w, res$w)
res$caste = 'None'
res$caste[res$padj < 0.05 &res$log2FoldChange < 0] = 'Worker'
res$caste[res$padj < 0.05 &res$log2FoldChange > 0] = 'Gyne'
table(res$caste)
boxplot(w ~ caste ,data = res)

boxplot(deg_overlap_w$w,deg_overlap_g$w , res$w)
wilcox.test(deg_overlap_g$w, deg_overlap_nde$w)
wilcox.test(deg_overlap_w$w, res$w)
wilcox.test(deg_overlap_g$w, deg_overlap_w$w)
cor.test(res$log2FoldChange ,res$w,method = 'k')
