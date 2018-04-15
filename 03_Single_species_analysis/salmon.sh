header=Aech #Only using A.echinatior as an example, we did the same for other species.
genome_index=ref_genome/${header}/${header} #Index genome using the output of Gemoma reannotation (and added ERCC if ERCC were added before cDNA library)
path_dir=raw_data/sorted/${header} #Directory where raw RNAseq files (Fastq) were stored.
for i in `ls $path_dir`;do
    if [ -d $path_dir/$i ];then
        echo "processing $i"
        nice salmon quant -i $genome_index -l A \
            -1 $path_dir/$i/*/*1.fq.gz \
            -2 $path_dir/$i/*/*2.fq.gz \
            --seqBias --gcBias \
            -p 10 -o quants/${i}
    fi
done
