library(tidyr)
library(ggplot2)
library(binom)
library(PropCIs)
gene_age_aech = read.table('Aech_age_odd.csv')
gene_age_aech$Freq = gene_age_aech$Freq+1
gene_age_aech_wide <- spread(gene_age_aech, Var1, Freq)

names(gene_age_aech_wide) = c('Age','Gyne','Major','Minor','NDE','Worker')
gene_age_aech_wide$Age = factor(gene_age_aech_wide$Age, levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
gene_age_aech_wide$Species = "A.echinatior"
gene_age_aech_wide[,c(8:10)] = binom.confint(gene_age_aech_wide$Gyne,rep(sum(gene_age_aech_wide$Gyne),8),conf.level = 0.9,methods = 'exact')[,c(4:6)]/
  rowSums(gene_age_aech_wide[,c(2:6)])*sum(gene_age_aech_wide[,c(2:6)])
gene_age_aech_wide[,c(11:13)] = binom.confint(gene_age_aech_wide$Major,rep(sum(gene_age_aech_wide$Major),8),conf.level = 0.9,methods = 'exact')[,c(4:6)]/
  rowSums(gene_age_aech_wide[,c(2:6)])*sum(gene_age_aech_wide[,c(2:6)])
gene_age_aech_wide[,c(14:16)] = binom.confint(gene_age_aech_wide$Minor,rep(sum(gene_age_aech_wide$Minor),8),conf.level = 0.9,methods = 'exact')[,c(4:6)]/
  rowSums(gene_age_aech_wide[,c(2:6)])*sum(gene_age_aech_wide[,c(2:6)])



gene_age_aech_wide_long = gather(gene_age_aech_wide,caste,portion_mean,mean,mean.1,mean.2,factor_key = T)
gene_age_aech_wide_long$portion_lower = gather(gene_age_aech_wide,caste,portion_lower,lower,lower.1,lower.2,factor_key = T)[,15]
gene_age_aech_wide_long$portion_upper = gather(gene_age_aech_wide,caste,portion_upper,upper,upper.1,upper.2,factor_key = T)[,15]
levels(gene_age_aech_wide_long$caste) = c('Gyne',"Major worker","Minor worker")

read_odd = function(input_file,species){
  gene_age_aech = read.table(input_file)
  gene_age_aech$Freq = gene_age_aech$Freq+1
  gene_age_aech_wide <- spread(gene_age_aech, Var1, Freq)
  names(gene_age_aech_wide) = c('Age','Gyne','NDE','Worker')
  gene_age_aech_wide$Age = factor(gene_age_aech_wide$Age, levels = c('neoptera','endopterygota','hymenoptera','aprocrita','aculeata','formicidae','myrmicinae','TRG'))
  gene_age_aech_wide$Species = species
  gene_age_aech_wide[,c(6:8)] = binom.confint(gene_age_aech_wide$Gyne,rep(sum(gene_age_aech_wide$Gyne),8),conf.level = 0.9,methods = 'exact')[,c(4:6)]/
    rowSums(gene_age_aech_wide[,c(2:4)])*sum(gene_age_aech_wide[,c(2:4)])
  gene_age_aech_wide[,c(9:11)] = binom.confint(gene_age_aech_wide$Worker,rep(sum(gene_age_aech_wide$Worker),8),conf.level = 0.9,methods = 'exact')[,c(4:6)]/
    rowSums(gene_age_aech_wide[,c(2:4)])*sum(gene_age_aech_wide[,c(2:4)])
  gene_age_aech_wide_long = gather(gene_age_aech_wide,caste,portion_mean,mean,mean.1,factor_key = T)
  gene_age_aech_wide_long$portion_lower = gather(gene_age_aech_wide,caste,portion_lower,lower,lower.1,factor_key = T)[,11]
  gene_age_aech_wide_long$portion_upper = gather(gene_age_aech_wide,caste,portion_upper,upper,upper.1,factor_key = T)[,11]
  levels(gene_age_aech_wide_long$caste) = c('Gyne',"Worker")
  return(gene_age_aech_wide_long[,c(1,5,10:13)])
}
Mpha_long = read_odd('Mpha_age_odd.csv','M.pharaonis')
Sinv_long = read_odd('Sinv_age_odd.csv','S.invicta')
Lhum_long = read_odd('Lhum_age_odd.csv','L.humile')
Lnig_long = read_odd('Lnig_age_odd.csv','L.niger')
Cbir_long = read_odd('../Cbir/Cbir_age_odd.csv','C.biroi')
Dqua_long = read_odd('../Dqua/Dqua_age_odd.csv','D.qua')



Aech_long = gene_age_aech_wide_long[,c(1,7,14:17)]
four_long = rbind(Mpha_long,Sinv_long,Lhum_long,Lnig_long,Aech_long)
four_long$Species = factor(four_long$Species,levels = c('A.echinatior','S.invicta','M.pharaonis','L.niger','L.humile'))
four_long$likelihood = log(four_long$portion_mean)
four_long$likelihood_lower = log(four_long$portion_lower)
four_long$likelihood_higher = log(four_long$portion_upper)
names(four_long)[3] = 'Caste-biased'
levels(four_long$Age) = c('(Pre-)Neoptera','Endopterygota','Hymenoptera','Aprocrita','Aculeata','Formicidae','Myrmicinae','Lineage specific')
four_long = four_long[-which(four_long$`Caste-biased` == 'Major worker'),]
four_long[four_long$`Caste-biased` == 'Minor worker','Caste-biased'] = "Worker"
four_long$likelihood_lower[four_long$likelihood_lower < -2] = -2
library(tidyverse)
a = ggplot(four_long,aes(x = Age, y = likelihood,group = `Caste-biased`))+
  geom_point()+
  geom_line(data=four_long,aes(colour = `Caste-biased`))+
  geom_ribbon(data=four_long,aes(ymin=likelihood_lower,ymax=likelihood_higher, fill = `Caste-biased`),alpha=0.2)+
  facet_wrap(~Species,ncol = 1,scales = 'free_y')+
 # ylim(c(-2,1.5))+
  geom_hline(yintercept = 0)+
  scale_fill_manual(values = c('red','blue'))+
  scale_color_manual(values = c('red','blue'))+
  theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1, size=10),legend.position="right",legend.direction = 'vertical')+
  ylab('Logarithmic Likelihood-Ratio') + xlab("Evolutionary Origin")
a
  
