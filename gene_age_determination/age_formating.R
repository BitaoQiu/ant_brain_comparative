library(readr)
#myproject <- read_delim("myproject.proteinortho","\t", escape_double = FALSE, trim_ws = TRUE)
myproject <- read_delim("myproject.poff","\t", escape_double = FALSE, trim_ws = TRUE)

names(myproject)
ant_ref = paste(c('Aech','Mpha','Sinv','Lnig','Lhum','Cbir','Dqua','Hsal_ant'),
                'fasta',sep = '.')
bee_ref = names(myproject)[grep('bee.fasta',names(myproject))]
wasp_ref = paste(c('Cflo','Csol','Nvit','Tpre'),'wasp.fasta',sep = '_')
hymenoptera_ref = "Oabi_wasp.fasta"
endopterygota_ref = paste(c('Dmel_fruitfly','Mcin_butterfly','Tcas_beetle'),
                          'fasta',sep ='.')
neoptera_ref = paste(c('Focc_thrip','Phum_louse','Znev_termite'),
                     'fasta',sep = '.')
sum(myproject[1,neoptera_ref] == '*') < length(myproject[1,neoptera_ref])
myproject[1,neoptera_ref] == '*'

ortholog_find = function(orth_data,lineage_input){
  return(apply(orth_data,1,FUN = function(x,lineage = lineage_input){
    return(sum(x[lineage] == '*') < length(x[lineage]))}))}
check_age = function(x){
  for(stage in seq(7,1,-1)){
    if (stage > 2){if (x[stage] == TRUE & sum(x[c(1,2)]) >= 1 & 
                       dbinom(sum(x[3:stage] == 1) ,size = length(x[3:stage]),prob = .9)> 0.01){
      return(names(x)[stage])
    }
    }
    if (stage == 2 & x[2] == 1){
      return(names(x)[stage])
    }
    if (stage == 1 & x[1] == 1){
      return(names(x)[stage])
    }
  }
  return('TRG')
}

check_age_bee = function(x){
  for(stage in seq(6,1,-1)){
    if (stage > 1){if (x[stage] == TRUE & sum(x[c(1)]) >= 1 & 
                       dbinom(sum(x[2:stage] == 1) ,size = length(x[2:stage]),prob = .9)> 0.01){
      return(names(x)[stage])
    }
    }
    if (stage == 1 & x[1] == 1){
      return(names(x)[stage])
    }
  }
  return('TRG')
}


head(ortholog_find(myproject,bee_ref))
####
Amel_data = myproject[-c(which(myproject$Amel_bee.fasta == '*'),
                         grep(',',myproject$Amel_bee.fasta)),]


Amel_age = data.frame(apoidea = ortholog_find(Amel_data,paste(c('Bter_bee','Hlab_bee','Mqua_bee'),'fasta',
                                                                 sep = '.')),
                      aculeata = ortholog_find(Amel_data,ant_ref),
                      aprocrita = ortholog_find(Amel_data,wasp_ref),
                      hymenoptera = ortholog_find(Amel_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Amel_data,endopterygota_ref),
                      neoptera = ortholog_find(Amel_data,neoptera_ref),row.names =  Amel_data$Amel_bee.fasta)
Amel_age = Amel_age*1
Amel_age_all = data.frame(age = apply(Amel_age,1,check_age_bee))
head(Amel_age_all)
table(Amel_age_all)
write.table(Amel_age_all,'Amel_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')

####
Aech_data = myproject[-c(which(myproject$Aech.fasta == '*'),
                        grep(',',myproject$Aech.fasta)),]


Aech_age = data.frame(myrmicinae = ortholog_find(Aech_data,paste(c('Mpha','Sinv'),'fasta',
                                           sep = '.')),
           formicidae = ortholog_find(Aech_data,c(ant_ref[c(4:8)])),
           aculeata = ortholog_find(Aech_data,bee_ref),
           aprocrita = ortholog_find(Aech_data,wasp_ref),
           hymenoptera = ortholog_find(Aech_data,hymenoptera_ref),
           endopterygota = ortholog_find(Aech_data,endopterygota_ref),
           neoptera = ortholog_find(Aech_data,neoptera_ref),row.names =  Aech_data$Aech.fasta)
Aech_age = Aech_age*1
tmp = sapply(Aech_data[,c(paste(c('Mpha','Sinv'),'fasta',sep = '.'),c(ant_ref[c(4:8)]), bee_ref,wasp_ref,hymenoptera_ref,
             endopterygota_ref, neoptera_ref)],FUN = function(x)
             {return(x == '*')})
plot(hclust(dist(t(as.matrix(tmp*1)))))
head(Aech_age)


Aech_age_all = data.frame(age = apply(Aech_age,1,check_age))
head(Aech_age_all)
table(Aech_age_all)
write.table(Aech_age_all,'Aech_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')

Mpha_data = myproject[-c(which(myproject$Mpha.fasta == '*'),
                         grep(',',myproject$Mpha.fasta)),]


Mpha_age = data.frame(myrmicinae = ortholog_find(Mpha_data,paste(c('Aech','Sinv'),'fasta',
                                                                 sep = '.')),
                      formicidae = ortholog_find(Mpha_data,c(ant_ref[c(4:8)])),
                      aculeata = ortholog_find(Mpha_data,bee_ref),
                      aprocrita = ortholog_find(Mpha_data,wasp_ref),
                      hymenoptera = ortholog_find(Mpha_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Mpha_data,endopterygota_ref),
                      neoptera = ortholog_find(Mpha_data,neoptera_ref),row.names =  Mpha_data$Mpha.fasta)
Mpha_age = Mpha_age*1
Mpha_age_all = data.frame(age = apply(Mpha_age,1,check_age))
table(Mpha_age_all)
write.table(Mpha_age_all,'Mpha_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')

#####
Sinv_data = myproject[-c(which(myproject$Sinv.fasta == '*'),
                         grep(',',myproject$Sinv.fasta)),]


Sinv_age = data.frame(myrmicinae = ortholog_find(Sinv_data,paste(c('Aech','Mpha'),'fasta',
                                                                 sep = '.')),
                      formicidae = ortholog_find(Sinv_data,c(ant_ref[c(4:8)])),
                      aculeata = ortholog_find(Sinv_data,bee_ref),
                      aprocrita = ortholog_find(Sinv_data,wasp_ref),
                      hymenoptera = ortholog_find(Sinv_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Sinv_data,endopterygota_ref),
                      neoptera = ortholog_find(Sinv_data,neoptera_ref),row.names =  Sinv_data$Sinv.fasta)
Sinv_age = Sinv_age*1
Sinv_age_all = data.frame(age = apply(Sinv_age,1,check_age))
table(Sinv_age_all)
write.table(Sinv_age_all,'Sinv_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')
#####
Lnig_data = myproject[-c(which(myproject$Lnig.fasta == '*'),
                         grep(',',myproject$Lnig.fasta)),]


Lnig_age = data.frame(myrmicinae = ortholog_find(Lnig_data,paste(c('Aech','Mpha','Sinv'),'fasta',
                                                                 sep = '.')),
                      formicidae = ortholog_find(Lnig_data,c(ant_ref[c(5:8)])),
                      aculeata = ortholog_find(Lnig_data,bee_ref),
                      aprocrita = ortholog_find(Lnig_data,wasp_ref),
                      hymenoptera = ortholog_find(Lnig_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Lnig_data,endopterygota_ref),
                      neoptera = ortholog_find(Lnig_data,neoptera_ref),row.names =  Lnig_data$Lnig.fasta)
Lnig_age = Lnig_age*1
Lnig_age_all = data.frame(age = apply(Lnig_age,1,check_age))
table(Lnig_age_all)
write.table(Lnig_age_all,'Lnig_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')

#####
Lhum_data = myproject[-c(which(myproject$Lhum.fasta == '*'),
                         grep(',',myproject$Lhum.fasta)),]


Lhum_age = data.frame(myrmicinae = ortholog_find(Lhum_data,paste(c('Aech','Mpha','Sinv','Lnig'),'fasta',
                                                                 sep = '.')),
                      formicidae = ortholog_find(Lhum_data,c(ant_ref[c(6:8)])),
                      aculeata = ortholog_find(Lhum_data,bee_ref),
                      aprocrita = ortholog_find(Lhum_data,wasp_ref),
                      hymenoptera = ortholog_find(Lhum_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Lhum_data,endopterygota_ref),
                      neoptera = ortholog_find(Lhum_data,neoptera_ref),row.names =  Lhum_data$Lhum.fasta)
Lhum_age = Lhum_age*1
Lhum_age_all = data.frame(age = apply(Lhum_age,1,check_age))
table(Lhum_age_all)
write.table(Lhum_age_all,'Lhum_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')
#####

Dqua_data = myproject[-c(which(myproject$Dqua.fasta == '*'),
                         grep(',',myproject$Dqua.fasta)),]


Dqua_age = data.frame(myrmicinae = ortholog_find(Dqua_data,c(ant_ref[c(8)])),
                      formicidae = ortholog_find(Dqua_data,paste(c('Aech','Mpha','Sinv','Lnig','Lhum','Cbir'),'fasta',
                                                                 sep = '.')),
                      aculeata = ortholog_find(Dqua_data,bee_ref),
                      aprocrita = ortholog_find(Dqua_data,wasp_ref),
                      hymenoptera = ortholog_find(Dqua_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Dqua_data,endopterygota_ref),
                      neoptera = ortholog_find(Dqua_data,neoptera_ref),row.names =  Dqua_data$Dqua.fasta)
Dqua_age = Dqua_age*1
Dqua_age_all = data.frame(age = apply(Dqua_age,1,check_age))
table(Dqua_age_all)
write.table(Dqua_age_all,'Dqua_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')
#####

Cbir_data = myproject[-c(which(myproject$Cbir.fasta == '*'),
                         grep(',',myproject$Cbir.fasta)),]


Cbir_age = data.frame(myrmicinae = ortholog_find(Cbir_data,paste(c('Aech','Mpha','Sinv','Lnig','Lhum'),'fasta',
                                                                 sep = '.')),
                      formicidae = ortholog_find(Cbir_data,c(ant_ref[c(7,8)])),
                      aculeata = ortholog_find(Cbir_data,bee_ref),
                      aprocrita = ortholog_find(Cbir_data,wasp_ref),
                      hymenoptera = ortholog_find(Cbir_data,hymenoptera_ref),
                      endopterygota = ortholog_find(Cbir_data,endopterygota_ref),
                      neoptera = ortholog_find(Cbir_data,neoptera_ref),row.names =  Cbir_data$Cbir.fasta)
Cbir_age = Cbir_age*1
Cbir_age_all = data.frame(age = apply(Cbir_age,1,check_age))
table(Cbir_age_all)
write.table(Cbir_age_all,'Cbir_age.txt',quote = F,row.names = T, col.names = T, sep = '\t')
