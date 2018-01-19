library(limma)
targets <- readTargets("E-MTAB-2621.sdrf.txt")
names(targets)[25] = 'FileName'

f <- function(x) as.numeric(x$Flags > -99)
RG <- read.maimages(targets, source="genepix",names = targets$Source.Name, green.only = T,wt.fun=f)

RG <- backgroundCorrect(RG, method="normexp", offset=50)
summary(RG$E)
names(RG$genes)
MA <- normalizeWithinArrays(RG, method="loess")
dim(RG$genes)
dim(RG$E)
RG$targets
