norm_reference = combat_edata_train - rowMeans(combat_edata_train)
norm_target = combat_edata_target - rowMeans(combat_edata_target)
share_gene = intersect(rownames(norm_reference),rownames(norm_target)) #Only select genes that shared between both reference and target species.
eigen_reference = princomp(norm_reference[share_gene,]) #Extracted the eigenvector from reference.

pc_data_target = data.frame(t(t(eigen_reference$scores[,c(1,2)]) %*% norm_target[share_gene,])/100) # Dimension reduction for target species based on the eigenvector of reference species' GRN (the 1st and 2nd PCs)
pc_data_target$species = filter_table$species
pc_data_target$caste = filter_table$caste
levels(pc_data_target$species) = c('A.echinatior','C.biroi','D.quadriceps','L.humile','L.niger',"M.pharaonis","S.invicta")
pc_data_target$species = factor(pc_data_target$species,levels = c('A.echinatior',"S.invicta","M.pharaonis",'L.niger','L.humile','C.biroi','D.quadriceps'))

levels(pc_data_target$caste)  = c("Gyne","Worker","Non-reproductive","Reproductive",'Worker') #Note we put minor worker and worker as the same caste, but it does not change the conclusion.
pc_data_target$caste = factor(pc_data_target$caste,levels = c("Gyne","Worker", "Reproductive","Non-reproductive"))

names(pc_data_target) = c('PC1','PC2','Species',"Caste")
ggplot(pc_data_target, aes(PC1, PC2, color=Caste, shape=Species,alpha = Species)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(15:17,0,1,2,5))+
  scale_alpha_manual(values= c(rep(0.2,5),1,1))+
  scale_color_manual(values = c('red','blue','purple','darkblue'))+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) 

