for f in *R1*.fastq*
do
i="${f%_*R1*.fastq*}"
cutadapt \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o "$i"_R1_trimmed.fastq.gz \
-p "$i"_R2_trimmed.fastq.gz \
"$i"_R1_001.fastq.gz "$i"_R2_001.fastq.gz
done
