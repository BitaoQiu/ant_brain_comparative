# Examine the influence of species identity and caste status on gene expression level.

# Step 1. Prepared for orthologous genes table across ant species (and included honeybee in later stage).
nice python3 format_gene_ortholog_table.py # Required ${species}.gff files and the output of Proteinortho.

# Step 2. Clustering and PCA analyses for all samples based on the orthologous genes' expression levels.
# See without_normalization.R for the clustering and PCA result. Note that we transformed the data first with log and quantile normalization.

# Step 3. Normalization for between species/colony variation, and the clustering and PCA results after normalization. Detailed explanation of method can see the method section.
# See with_normalization.R for the clusterin and PCA results after normalization of species/colony variation.