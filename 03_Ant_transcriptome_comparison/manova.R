library("dplyr")
library(car)
library(nnet)
head(exp)
head(sampleTable)
sampleTable$lineage = factor(c(rep('m',8),rep('l',8),rep('m',18)))
sampleTable$worker = factor(c(rep('p',8),rep('m',18),rep('p',8)))

sampleTable$habit = factor(c(rep('a',12),rep('l',26),rep('o',20)))

plot(density(exp[,10]))


exp2 = exp[which(rownames(exp) %in% deg_union),]
exp2 = combat_edata
svd1 = svd(exp2 - rowMeans(exp2))
plot(svd1$d^2/sum(svd1$d^2), xlim = c(0, 15), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")

plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     pch = as.numeric(sampleTable$species)+13, col = c('red','blue'))
pc_data = svd1$v
colnames(pc_data) = paste('pc',c(1:34),sep = '')
rownames(pc_data) = rownames(sampleTable)
pc_data = as.data.frame(cbind(pc_data,sampleTable ))
Manova(lm1 <- lm(cbind(pc1,pc2) ~  caste*species , data=pc_data), type="III")

summary(aov(cbind(pc1,pc2) ~ caste*species, data = pc_data))
summary(aov(cbind(pc1,pc2) ~ caste*lineage, data = pc_data))

Anova(lm1 <- glm( pc2 ~  caste*lineage, data=pc_data) , type="III")

summary(res.man)
summary.aov(res.man)

lm1 <- lm(pc1 ~  caste*species , data=pc_data)
lm2 <- lm(pc1 ~  caste*lineage , data=pc_data)
BIC(lm1,lm2)
