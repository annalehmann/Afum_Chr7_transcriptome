# build index for A1163
hisat2-build -p 16 FungiDB-68_AfumigatusA1163_Genome.fasta AfumigatusA1163

for f in *R1_trimmed.fastq*
do
i="${f%_*R1_trimmed.fastq*}"
hisat2 -q -x AfumigatusA1163 \
-1 "$i"_R1_trimmed.fastq* \
-2 "$i"_R2_trimmed.fastq* | \
samtools sort -o "$i"_HISAT2_aln.bam
done
