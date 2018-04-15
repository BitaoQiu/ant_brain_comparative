#Using the gene structures of the reference species, we re-annotated the seven ant species + other insect species (For later gene age calculation).
gemoma_jar=GeMoMa-1.4-beta2.jar
bam_dir=#directory that stored the RNAseq bam files, aligned with Hisat2 (whole genome alignment)

#Step 1. For the seven ant species, for each species, we first extracted the exon information based on RNAseq data.
bam_list=`ls $bam_dir/*/*sorted.bam | paste -d"," -s| sed "s/\,/\ m=/g"`
nice java -jar $gemoma_jar CLI ERE m=$bam_list c=true

#Step 2. We then independently reannotated the target genome with the reference genomes and integrated the RNAseq evidence:
genome=#The NCBI genome sequence of target species
tmp_ref=#The prepared reference genome from 1.1.1.Prepare_reference_genomes.sh
intron_gff=#Output from Step 1.
coverage_file=#Output from Step 1.
for out in `ls $tmp_ref`;do #For each reference genome
    if [ -d $tmp_ref/$out ];then
        mkdir $out
        makeblastdb -out ${out}/blastdb -hash_index -in ${genome} -title "target" -dbtype nucl #Prepare for blast
        for i in `ls $tmp_ref/$out/split-*.fasta`;do #split-*fasta correponding the split fasta files from preparation, so that can perform blast in parallel.
            echo $i
            out_1=`basename ${i/split/tblasn}`
            out_file=$out/${out_1/fasta/txt}
            nice tblastn -query $i -db ${out}/blastdb -evalue 100.0 -out $out_file -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1 && java -jar ${gemoma_jar} CLI GeMoMa t=$out_file c=$i i=$intron_gff a=$tmp_ref/$out/assignment.tabular coverage='UNSTRANDED' coverage_unstranded=$coverage_file tg=${genome} outdir=${out} ct=0.4 p=10 # This step can run in parallel.
        done
    fi
done

#Step 3. Concatenated the annotations of the split fasta files (output of step2) into one.
for i in `ls`;do
    if [ -d $i ];then
        cat $i/predicted_annotation.gff >> $i/predicted_annotation_all.gff
        cat $i/predicted_protein.fasta >> $i/predicted_protein_all.fasta
        for j in {1..9};do
            cat $i/predicted_annotation_$j.gff >> $i/predicted_annotation_all.gff
            cat $i/predicted_protein_$j.fasta >> $i/predicted_protein_all.fasta
        done
    fi
done

#Step 4. Because we reannotated the target genome indpendently from multiple reference genomes, to trace the origin of gene model, we re-named the gff file by adding the reference species as prefix in GeneID.
for i in `ls`;do
    if [ -d $i ];then
        sed "s/;ref-gene=/;ref-gene=${i}_/g" ${i}/predicted_annotation_all.gff | sed "s/ID=/ID=${i}_/g" | sed "s/;Parent=/;Parent=${i}_/g" > ${i}/predicted_annotation_all_named.gff
        sed "s/>/>${i}_/g" ${i}/predicted_protein_all.fasta > ${i}/predicted_protein_all_named.fasta
    fi
done

#Step 5. Combined the annotations and perform filtering (e.g. Overlap, missing stop codon or start codon, and RNAseq evidence), see http://www.jstacs.de/index.php/GeMoMa for details.
gff_file=`ls */predicted_annotation_all_named.gff | paste -d"," -s| sed "s/\,/\ g=/g"`
nice java -jar $gemoma CLI GAF g=$gff_file e=0.1
sed 's/prediction/mRNA/g' filtered_predictions.gff > filtered_predictions_named.gff

# The filtered_predictions_named.gff can be used from transcriptome alignment or orthologous gene detection.

# For other insect species, we performed the re-annotation for each species using the same steps as above, except without using the RNA data.

