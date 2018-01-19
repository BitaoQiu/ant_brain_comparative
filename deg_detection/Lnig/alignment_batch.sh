header=Lnig
genome_index=/isdata/cse/bitao/gene_age_mar_2017/ref_genome/${header}/${header}
path_dir=/isdata/cse/bitao/gene_age_mar_2017/raw_data/${header}
for i in `ls $path_dir`;do
    if [ -d $path_dir/$i ];then
        echo "processing $i"
        nice salmon quant -i $genome_index -l A \
            -1 $path_dir/$i/*_1.fastq \
            -2 $path_dir/$i/*_2.fastq \
            --seqBias --gcBias \
            -p 10 -o quants/${i}
    fi
done
