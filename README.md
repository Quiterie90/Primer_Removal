# Primer_Removal
## Input files
This script takes as input two files:
- a fastq file containing joined paired-ends reads
- a primers file
The primers file consists in a two-lines file where the first line corresponds to the sequence of the forward primer in 5'-> 3' orientation and the second line correcponds to the sequence of the reverse primer in 5'-> 3' orientation.
This file is generated by the user. 

## How to run the script
python ./remove_primers_fastq.py FileName.fastq PrimerName.primers

## How does it work
The script retains the reads that present the exact sequence of the forward primer at the beginning of the read and the exact sequence of the reverse primer at the end of the read and trims both extremities.

## Output files
Two fastq output files are generated:
- one containing the trimmed reads -> FileName_trimmed.fastq
- one containing the discarded reads -> FileName_failed.fastq
