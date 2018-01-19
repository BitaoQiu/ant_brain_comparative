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
t2gene[grep('odorant',t2gene$blast),]
head(t2gene)
deg_overlap_w = round(log2_deg[which(rowSums(padj_deg[,c(2:5)] < .1) == 4 & rowSums(log2_deg[,c(2:5)] < 0) == 4),c(1:5)],2)
deg_overlap_w$dmel = GO_annotation$name[match(rownames(deg_overlap_w),GO_annotation$V2)]
deg_overlap_w$id = t2gene$blast[match(rownames(deg_overlap_w),t2gene$V2)]
deg_overlap_w
deg_overlap_g = round(log2_deg[which(rowSums(padj_deg[,c(2:5)] < .1) >= 4 & rowSums(log2_deg[,c(2:5)] > 0) == 4),c(1:5)],2)
deg_overlap_g$dmel = GO_annotation$name[match(rownames(deg_overlap_g),GO_annotation$V2)]
deg_overlap_g$id = t2gene$blast[match(rownames(deg_overlap_g),t2gene$V2)]
deg_overlap_g[,c(-1)]

geneid = 'Aech_gene_5668'
check_exp(ortholog_exp,geneid,title = t2gene[t2gene$V2 == geneid,'blast'])#+ geom_abline(slope = sqrt(2),intercept = 0)+ geom_abline(slope = sqrt(.5),intercept = 0)

check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_13482','Aech_gene_15227','Aech_gene_9986'),'Genes with worker biased expression',
               gene_name = c('Calcium/calmodulin-dependent protein kinase',
                             'Ras-related and estrogen-regulated growth inhibitor-like protein',
                             'gamma-aminobutyric acid receptor subunit'))
check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_13364','Aech_gene_5319','Aech_gene_14561'),'Genes with gyne biased expression',
               gene_name = c('neuroparsin-A-like isoform',
                             'insulin-like growth factor-binding protein',
                             'arrestin homolog'))

check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_13482','Aech_gene_9986'),'Behaviours related',
               gene_name = c('Calcium/calmodulin-dependent protein kinase',
                             'gamma-aminobutyric acid receptor subunit'))
check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_13364','Aech_gene_5533','Aech_gene_5319','Aech_gene_15120'),'Juvenile hormone related',
               gene_name = c('neuroparsin-A-like isoform',
                             'vitellogenin-3-like','insulin-like growth factor-binding protein','insulin-like growth factor II'))
check_exp_list(ortholog_exp,gene_id_list = c('Aech_gene_15227','Aech_gene_5319'),'Juvenile hormone related',
               gene_name = c('Ras-related and estrogen-regulated growth inhibitor-like protein',
                             'insulin-like growth factor-binding protein'))

t2gene[which(t2gene$V1 == 'AECH_RNA9173_R0'),]
log2_deg[which(rownames(log2_deg) == geneid),]
padj_deg[which(rownames(padj_deg) == geneid),]
