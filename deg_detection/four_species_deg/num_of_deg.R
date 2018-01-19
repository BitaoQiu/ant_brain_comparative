#Check DEGs
cut_off = 0.2
deg_overlap_w = round(log2_deg[which(rowSums(padj_deg[,c(2:5)] < .05) >= 2 & rowSums(log2_deg[,c(2:5)] < -cut_off) == 4),c(1:5)],4)

deg_overlap_w$dmel = GO_annotation$name[match(rownames(deg_overlap_w),GO_annotation$V2)]
deg_overlap_w
deg_overlap_g = round(log2_deg[which(rowSums(padj_deg[,c(2:5)] < .05) >= 2 & rowSums(log2_deg[,c(2:5)] > cut_off) == 4),c(1:5)],4)
deg_overlap_g$dmel = GO_annotation$name[match(rownames(deg_overlap_g),GO_annotation$V2)]

dim(rbind(deg_overlap_w,deg_overlap_g))
deg_overlap_g

dim(log2_deg[which(log2_deg$Lhum_log2FoldChange < -0.6 & padj_deg$L.humile < 0.05),])
head(res)
