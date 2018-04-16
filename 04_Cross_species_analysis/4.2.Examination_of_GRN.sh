# The existence of reference genetic regulatory network (GRN) in target species can be examined by first extracting the eigenvectors (weighted combination of expression levels of genes associated with GRN) from reference species using SVD; then projecting the expression levels of target species onto the eigenvectors. In this sense, we are examining the expression level of GRN related genes in target species, if GRN exists in target species, it should exhibit similar expression pattern with the reference species.

# Step 1: Normalized the transcriptome of reference species (five ants with permanent castes) and target species (queenless ants, honey bee, or poliste wasp), independently.
# See salmon_reference.R and salmon_target.R. In the salmon_target.R, we are only using the queenless ants as example, but the same principles were used in both the honey bee and poliste wasp.

# Step 2: Extracted the eigenvectors from normalized transcriptome of reference species (representing GRN) using SVD and projected normalized transcriptome of target species onto the eigenvectors.
# See test_GRN.R

