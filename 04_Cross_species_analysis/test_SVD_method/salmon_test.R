tmp = combat_edata_train - rowMeans(combat_edata_train)
tmp_2 = combat_edata - rowMeans(combat_edata)
share_gene = intersect(rownames(tmp),rownames(tmp_2))
tmp_1 = princomp(tmp[share_gene,])

pc_data_test = data.frame(t(t(tmp_1$scores[,c(1,2)]) %*% tmp_2[share_gene,])/100)
#pc_data_test = data.frame(t(t(abs(tmp_1$scores[,c(1,2)])) %*% abs(tmp_2[share_gene,])/100))
pc_data_test$species = filter_table$species
pc_data_test$caste = filter_table$caste
levels(pc_data_test$species) = c('A.echinatior','C.biroi','D.quadriceps','L.humile','L.niger',"M.pharaonis","S.invicta")
pc_data_test$species = factor(pc_data_test$species,levels = c('A.echinatior',"S.invicta","M.pharaonis",'L.niger','L.humile','C.biroi','D.quadriceps'))

levels(pc_data_test$caste)  = c("Gyne","Worker","Non-reproductive","Reproductive",'Worker')
pc_data_test$caste = factor(pc_data_test$caste,levels = c("Gyne",'Worker', "Reproductive","Non-reproductive"))
plot_alpha = rep(0.3,7)
#plot_alpha[c(2,3)] = 1
plot_alpha[n] = 1

#percentVar <- round(100 * attr(pc_data_test, "percentVar"))
names(pc_data_test) = c('PC1','PC2','Species',"Caste")
ggplot(pc_data_test, aes(PC1, PC2, color=Caste, shape=Species,alpha = Species)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15:17,0,1,2,5))+
  scale_alpha_manual(values=plot_alpha)+
  scale_color_manual(values = c('red','blue','purple','black'))+
  #scale_x_continuous(trans = 'reverse')+
  xlab(paste0("PC1 (",percentVar_train[1],"%)")) +
  ylab(paste0("PC2 (",percentVar_train[2],"%)"))+
  labs(title = levels(pc_data_test$Species)[n])+
  coord_fixed()

