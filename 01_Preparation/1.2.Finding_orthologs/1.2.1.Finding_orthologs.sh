#Step 1. For each target genome, extracted the longest isoform for each gene, and prepared the synteny information, based on the re-annotation output.
ref_dir=#Output of 1.1.2, for all the species.
for i in `ls $ref_dir`;do #Do it for all the target species, $i is the corresponding species directory
    if [ -d $ref_dir/$i ];then
        python3 find_longestisoform_gemoma.py ${ref_dir}/${i}/filtered_predictions_named.gff ${ref_dir}/${i}/proteins.fasta $i
    fi
done

#Step 2. Detected orthologs using proteinortho (See https://www.bioinf.uni-leipzig.de/Software/proteinortho/)
nice proteinortho5.pl -synteny -singles -verbose *fasta

