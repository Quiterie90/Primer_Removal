# Primer_Removal
This script takes as input two files:
- a fastq file containing joined paired-ends reads
- a primers file
The primers file consists in a two-lines file where the first line corresponds to the sequence of the forward primer in 5'-> 3' orientation and the second line correcponds to the sequence of the reverse primer in 5'-> 3' orientation.

The script retains only the reads that present the exact sequence of the forward primer at the beginning of the read and the exact sequence of the reverse primer at the end of the read.

Trimmed reads are written in a new fastq file (FileName_trimmed.fastq) and the others are written in an other fastq file (FileName_failed).
