genome=ncbi_aech_withERCC_hisat/ncbi_aech_withERCC #Genome of A.echinatior (from NCBI) + ERCC sequences), indexed for Hisat2
path_dir=raw_data/sorted/Aech #Raw data of A.echinatior (fastq)
for i in `ls $path_dir`;do
    if [ -d $path_dir/$i ];then
        file_1=`ls $path_dir/$i/*/*_1.fq.gz | paste -d',' -s`
        file_2=`ls $path_dir/$i/*/*_2.fq.gz | paste -d',' -s`
        out_put_dir=`basename $i`
        mkdir $out_put_dir
        nice hisat2 -x $genome -1 $file_1 -2 $file_2 -S $out_put_dir/$out_put_dir.sam -p 10 --dta 2>$out_put_dir/mapping.txt
    fi
done


