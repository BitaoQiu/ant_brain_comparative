1) format_gene_ortholog_table.ipynb will format the input from Proteinortho, and generate GeneID orthologous table, as input for age_formating.R.
2) age_formating.R takes the input from format_gene_ortholog_table.ipynb, and determined the age of genes based on species phylogenetic relations, loss rate can be adjusted with "check_age" function, which can be used as input for /deg_detection/Aech/gene_network.R and further generated input for gene_age.R
3) For gene_age.R, we examined the different loss rate (5%, 10%, and 20%), and checked the association of gene age and log-likelihood ratio to be recruited as DEGs.

