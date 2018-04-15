# Using the output of Proteinortho (existence of orthologs across insects), we first formated it into ortholog matrix:
# Step 1. Formatted the output of Proteinortho, and reported the occurance of target gene orthologs across lineages (neoptera, endopterygota, hymenoptera, aprocrita, aculeata, formicidae and myrmicinae):
nice R CMD BATCH age_formatting.R # This will output the formatted gene age table for each target species.

# Step 2. Combined gene age and caste-biased expression information:
nice R CMD BATCH gene_age_summary.R #Only showing one species (M.pharaonis) to illustrate.

# Step 3. Visualization of the relationship between gene age and caste-biased expression:
nice R CMD BATCH gene_age_plot.R



