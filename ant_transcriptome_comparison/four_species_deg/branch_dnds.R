w = read.table('b_free_dnds.txt')
res$w_0 = w$V2[match(rownames(res),w$V1)]
res$w_1 = w$V3[match(rownames(res),w$V1)]
dnds = res[!is.na(res$w_0),]
w[w$V1 == 'Aech_gene_7096',]

deg_overlap_w$w_0 = w$V2[match(rownames(deg_overlap_w),w$V1)]
deg_overlap_w$w_1 = w$V3[match(rownames(deg_overlap_w),w$V1)]
deg_overlap_g$w_0 = w$V2[match(rownames(deg_overlap_g),w$V1)]
deg_overlap_g$w_1 = w$V3[match(rownames(deg_overlap_g),w$V1)]
boxplot(deg_overlap_w$w_0,deg_overlap_w$w_1,deg_overlap_g$w_0 , deg_overlap_g$w_1, dnds$w_0,dnds$w_1, outline = F,
        names= c(paste(rep(c('Dmel','ant'),3),rep(c('Worker','Gyne','All'),each = 2),sep ='-')), main = 'Branch model (ant-fruitfly)')
#boxplot(deg_overlap_w$w_0,deg_overlap_g$w_0 , dnds$w_0,outline = F)
res$caste = 'None'
res$caste[res$padj < 0.05 &res$log2FoldChange < 0] = 'Worker'
res$caste[res$padj < 0.05 &res$log2FoldChange > 0] = 'Gyne'
table(res$caste)
boxplot(dnds$caste == w_0,dnds$w_1,ylim = c(0,1))

boxplot(deg_overlap_w$w,deg_overlap_g$w , res$w)
wilcox.test(deg_overlap_g$w, res$w)
wilcox.test(deg_overlap_w$w, res$w)
wilcox.test(deg_overlap_g$w_1, deg_overlap_w$w_1)
cor.test(res$log2FoldChange ,res$w,method = 'k')
