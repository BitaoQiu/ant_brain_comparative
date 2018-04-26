# Codes for comparative ant brain transcriptomes
## 01. Preparation
# In this folder, we provided scripts for re-annotating genomes (1.1) with GeMoMwa and finding orthologs across seven ant species (1.2) with RBH and synteny information. The re-annonated genomes are used for downstream analayses, including gene age determination and expression quantification, and the ortholog information is used for cross species comparison.

## 02. Gene age determination
# We provided R and python scripts that used for gene age determination and visualization.

## 03. Single species analysis
# We provided bash and R script for gene expression quantification (Hisat+Stringtie and Salmon), sample quality determination, and caste differentially expressed genes (DEGs) identification.

## 04. Cross species comparison
# We provided R scripts for comparing orthologous gene expression across species (4.1), normalizing for species/colony identity on gene expression (4.1), using SVD to examine the genetic regulatory network in queenless ants, honeybee and paper wasp (4.2), constructing cross-species co-expression network (4.3), and identifying DEGs across ant species (4.4). 
_Note:_ *04_Cross_species_analysis/test_SVD_method* directory includes codes for testing SVD method (using leave-one-out Jackknife).

