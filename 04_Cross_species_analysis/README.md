# Guideline:
1. To examine the orthologous gene expression level across ant species, run: *without_normalization.R* for the clustering and PCA result.

2. To examine the expression level after normalizing for species/colony variation, run: *with_normalization.R*, note that it uses the output from *without_normalization.R*

3. To examine the existence of GRN in the two queenless ant species, 1) run *salmon_reference.R* to get normalized gene expression data of five typical ant species, 2) run *salmon_target.R* to get normalized gene expression data of the two queenless ant species, and 3) run *test_GRN.R* to extract the eigenvector (representing GRN) from the typical ant species using **SVD** and apply it to the queenless ant species. This will examine the GRN in the queenless ant species. 

4. To construct cross species co-expression network, run *network_five_ants.R*. All expression levels were normalized with colony identity.

5. To identify genes differentially expressed between castes across five ants species, run *cross_species_deg.R*

*Note*: To validate the **SVD** method, go to [test_SVD_method](https://github.com/StanQiu/ant_brain_comparative/tree/master/04_Cross_species_analysis/test_SVD_method), then run *salmon_trainning.R*, *salmon_all.R* and *salmon_test.R*.

# Input:
Input files for these analyses are available in [https://sid.erda.dk/sharelink/HfHUjTGrnv](https://sid.erda.dk/sharelink/HfHUjTGrnv).

# Requirements:
* R version: > 3
* R packages: devtools, Biobase, preprocessCore, tximport, DESeq2, RColorBrewer, pheatmap, ggplot2, sva, readr, rjson
