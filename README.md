# Codes for comparative ant brain transcriptomes
## Gene_age_determination
Used it to constructed gene age based on their ortholog occurrence on insect phylogentic tree.

1) *format_gene_ortholog_table.ipynb* will format the input from Proteinortho, and generate GeneID orthologous table, as input for age_formating.R.

2) *age_formating.R* takes the input from *format_gene_ortholog_table.ipynb*, and determined the age of genes based on species phylogenetic relations, loss rate can be adjusted with "check_age" function, which can be used as input for */Deg_detection/SPECIES/gene_network.R* and further generated input for *gene_age.R*

3) For *gene_age.R*, we examined the different loss rate (5%, 10%, and 20%), and checked the association of gene age and log-likelihood ratio to be recruited as DEGs.

## DEG_detection
Used it to detect differentially expressed genes between castes/phenotypes, separated for each species.

1) *find\_deg.R*/*salmon.R* used the output from _Salmon_ to detect DEGs with _DESeq2_.

2) *gene_network.R* constructed DEGs and gene age table (based on orthologous information).

## Ant_transcriptome_comparison
Used it to detect conserved caste differentially expressed genes across ant species, adjust for species/colony identity for gene expression, and test conservation of ant caste genetic regulatory network in queenless ants, honey bee, and social wasp.

1) *deg_full_model.R* will detect differentially expressed genes between castes across five ant species, using GLM model from *DESeq2*.

2) *salmon_all.R* will adjust for the species/colony identity for gene expression (caste specific gene expression), and visualize with *ggplot2* package in R.

3) *manova.R* will test the contribution of species/sub-family identity, caste status, and their interaction to the PCA data, using Multivariate analysis of variance.

4) *salmon_trainning.R* and *salmon_test.R* will test the existence of caste GRN in other species: train with *salmon_trainning.R* (using five ant species with normal caste differentiation as training data), then use *salmon_all.R* to include target species and test with *salmon_test.R*, which will use SVD to extract eigenvectors from training data and project data from target species onto the trained eigenvectors.

_Note:_ *SVD_test* directory includes codes for testing SVD method (using leave one out Jack knife).

