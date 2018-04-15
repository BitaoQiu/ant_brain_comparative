# We first reannotated the genomes of the seven ant species with the NCBI annotation of these seven species + honeybee +  Nasonia vitripennis + Drosophila.
#Step1. Extract the gene structure information from each species based on corresponding NCBI gff file:
annotation=#NCBI gff file
reference=#Genome sequence from NCBI.
out=#output directory
gemoma_jar=GeMoMa-1.4-beta2.jar #Find it in http://www.jstacs.de/index.php/GeMoMa
nice java -jar $gemoma_jar CLI Extractor a=$annotation g=$reference outdir=$out
#Step1.2. Split the fasta file into 10 files, in order to run Blast in parallel. (Optional)
nice java -cp $gemoma_jar projects.FastaSplitter cds-parts.fasta 10 "_"

