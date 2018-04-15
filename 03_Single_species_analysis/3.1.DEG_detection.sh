# Performed transcriptome alignment with re-annotated genomes (+ ERCC for A.echinatior, M.pharaonis and L.humile), and filtered samples based on ERCC and within species Spearman correlation coefficient. Then detected differentially expressed genes (DEGs) between castes for each species.

# Step 1.1 Transcriptome alignment with re-annotated genome (+ERCC) using Salmon (see https://combine-lab.github.io/salmon/):
# See salmon.sh for A.echinatior as an example.

# Step 1.2. Expression matrix using Hisat and StringTie (See https://www.nature.com/articles/nprot.2016.095):
# 1.2.1: Alignment with Hisat2:
# See hisat.sh for A.echinatior as an example.
# 1.2.2: Expression level quantification with Stringtie:
# See stringtie.sh for A.echinatior as an example.

# Step 2. Filtered samples based on ERCC or within species correlation coefficient:
# See check_quality.R for A.echinatior as an example.

# Step 3. Detected caste differentially expressed genes (DEGs) using DESeq2
# See find_deg.R for A.echinatior as an example.
