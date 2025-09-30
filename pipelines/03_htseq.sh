for file in *_HISAT2_aln.bam
do
base_name=${file%_HISAT2_aln.bam}
htseq-count \
-t protein_coding_gene \
-i ID \
-f bam \
$file \
FungiDB-68_AfumigatusA1163.gff > \
"$base_name"_counts_htseq.tsv
done
