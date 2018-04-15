sp=Aech
bam_dir=alignment/hisat_withERCC/$sp #Output of Hisat2
gff=ref_genome/$sp/filtered_predictions_named.gff #Output of Gemoma
gff2=ERCC/ERCC.gff #Gff for ERCC.fa
cat $gff $gff2 > stringtie_withERCC/$sp/gff_combined.gff #New gff for stringtie
gff_combined=stringtie_withERCC/$sp/gff_combined.gff
for i in `ls $bam_dir`;do
    if [ -d $bam_dir/$i ];then
        mkdir $i
        cd $i
        nice stringtie -e -B -o out.gtf -p 10 -G $gff_combined -A gene_abund.tab -C cov_refs.gtf $bam_dir/$i/${i}.sorted.bam 
        cd ..
    fi
done
