library(tximport)
library('DESeq2')
mpha_files = c(paste('Mpha',c(1:6,'7x',8,'9x',10),sep ='_'))#, paste('Mpha',c('7x','9x'),sep ='_'))

files <- file.path('quants', mpha_files, "quant.sf")
names(files) <- mpha_files
tx2gene <- read.csv('Mpha_gemoma_t2g.txt',header = F, sep = '\t',stringsAsFactors = F)
tx2gene$V1 = toupper(tx2gene$V1)
files
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
abundance = txi.salmon$abundance

# Between gyne and workers:
batch = factor(c(rep(c('gyne','worker'),5)))#,rep('gyne',2)))
colony = factor(rep(c(1:5),each = 2))
sampleTable <- data.frame(condition = batch, colony = colony)#, type = factor(c(rep('A',10),rep('B',2))))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition + colony)
#dds <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition','gyne','worker'),alpha = 0.05)
summary(res)
res = lfcShrink(dds, contrast=c("condition","gyne","worker"), res=res)
write.table(res, "Mpha_gyne_worker.csv",quote = F, row.names = T, col.names = T, sep = '\t')

######

