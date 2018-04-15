library(base)
library(tximport)
library('DESeq2')
aech_files = c(paste('Aech',c(4:9,'10x',11,'12x','13x',14,15),sep ='_')) #Selected samples after filtering.
files <- file.path('quants', aech_files, "quant.sf") #Output of salmon
names(files) <- aech_files
tx2gene <- read.csv('Aech_gemoma_t2g.txt',header = F, sep = '\t',stringsAsFactors = F) #Transcript and gene ID relationship, can be self generated from the output of Gemoma
tx2gene$V1 = toupper(tx2gene$V1)
#########
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,
                       countsFromAbundance = 'lengthScaledTPM')
cor(txi.salmon$abundance,method = 'p')[seq(3,12,3),seq(3,12,3)] #Examine the cross sample correlation coefficient.
batch = factor(c(rep(c('gyne','major','minor'),4)))
colony = factor(c(rep(c(1:4),each = 3)))
sampleTable <- data.frame(condition = batch, colony = colony)
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, 
                                ~condition  + colony) #Model the gene expression level with the influence of caste and colony.
dds <- DESeq(dds)
res_gb <- results(dds, contrast=c("condition","gyne","major"),alpha = 0.05) #Note: For Honeybee, using: alpha = 1e-3 and lfcThreshold = log2(1.5) as strict criteria of DEG detection, see Supplementary Text for explanation.
res_gs <- results(dds, contrast=c("condition","gyne","minor"),alpha = 0.05) #Detect DEGs between gyne and small workers.

