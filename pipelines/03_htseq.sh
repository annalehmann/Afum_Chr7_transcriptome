for f in *_HISAT2_aln.bam
do
i=${f%_HISAT2_aln.bam}
htseq-count \
-t protein_coding_gene \
-i ID \
-f bam \
$f \
FungiDB-68_AfumigatusA1163.gff > \
"$i"_counts_htseq.tsv
done
