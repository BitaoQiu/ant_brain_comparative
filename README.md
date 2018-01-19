# ant_brain_comparative
# Codes for comparative ant brain transcriptomes
## Gene_age_determination
Used it to constructed gene age based on their ortholog occurrence on insect phylogentic tree.

1) *format_gene_ortholog_table.ipynb* will format the input from Proteinortho, and generate GeneID orthologous table, as input for age_formating.R.

2) *age_formating.R* takes the input from *format_gene_ortholog_table.ipynb*, and determined the age of genes based on species phylogenetic relations, loss rate can be adjusted with "check_age" function, which can be used as input for */Deg_detection/SPECIES/gene_network.R* and further generated input for *gene_age.R*

3) For *gene_age.R*, we examined the different loss rate (5%, 10%, and 20%), and checked the association of gene age and log-likelihood ratio to be recruited as DEGs.

## DEG_detection
Used it to detect differentially expressed genes between castes/phenotypes, separated for each species.

*find\_deg.R*/*salmon.R* used the output from _Salmon_ to detect DEGs with _DESeq2_.

*gene_network.R* constructed DEGs and gene age table (based on orthologous information).

## Ant_transcriptome_comparison
