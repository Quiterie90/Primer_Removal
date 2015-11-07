#!/usr/bin/env python

import sys 

# Reading files
fastq_file = open(sys.argv[1], "r")
fastq_file_lines = fastq_file.read().split("\n")
fastq_file.close()

# Primers
primers_file = open(sys.argv[2], "r")
primers = primers_file.read().split("\n")
f_primer = primers[0]
r_primer = primers[1]
primers_file.close()

out_file_name = sys.argv[1].split('.fastq')[0] + "_trimmed.fastq"
out_file_name2 = sys.argv[1].split('.fastq')[0] + "_failed.fastq"

# This function allows the conversion of a DNA sequence 
# (in 5'->3') in its reverse complement (in 3'-> 5')
def do_complementary_r_primer(sequence):
	rc_primer = []
	for nucleotide in sequence:
		if nucleotide == 'A':
			rc_primer.append('T')
		elif nucleotide == 'T':
			rc_primer.append('A')
		elif nucleotide == 'C':
			rc_primer.append('G')
		elif nucleotide == 'G':
			rc_primer.append('C')
		elif nucleotide == 'W':
			rc_primer.append('W')
		elif nucleotide == 'S':
			rc_primer.append('S')
		elif nucleotide == 'Y':
			rc_primer.append('R')
		elif nucleotide == 'R':
			rc_primer.append('Y')
		else:
			print 'error, this is not a nucleotide!'

	return rc_primer[::-1]

# Comparison between primer and the fasta sequence at each position 
# Takes into account the ambiguous bases
def nucleotide_comparison(nucleotide_primer, nucleotide_sequence):
	if nucleotide_primer == nucleotide_sequence:
		return True
	else:
		if nucleotide_primer == 'W':
			if nucleotide_sequence in ('A','T'):
				return True
		elif nucleotide_primer == 'S':
			if nucleotide_sequence in ('G','C'):
				return True
		elif nucleotide_primer == 'Y':
			if nucleotide_sequence in ('C','T'):
				return True
		elif nucleotide_primer == 'R':
			if nucleotide_sequence in ('A','G'):
				return True
		elif nucleotide_sequence in ('A','C','T','G','N'):
			return False
		else:
			print "Error, problem with the nucleotide comparison" + nucleotide_sequence
	return False

# Main
rc_primer = do_complementary_r_primer(r_primer)

# Open output file
output_file = open(out_file_name, "w+")
output_file2 = open(out_file_name2, "w+")

# Length
seq_length_total = 0
cpt_sequences_matched = 0
cpt_sequences_unmatched =0
cpt_bad_f_primers = 0
cpt_bad_r_primers = 0
for i in xrange(0,len(fastq_file_lines) - 1,4):
	sequence = fastq_file_lines[i+1]
	result = True
	# print sequence

	# Checking forward primer
	for j in xrange(len(f_primer)):
		if not nucleotide_comparison(f_primer[j], sequence[j]):
			cpt_bad_f_primers += 1
			result = False
			break
	
	# If matching the primer, it is not necessary to check the reverse primer
	# Checking the reverse complementary primer
	if result:
		for j in xrange(len(rc_primer)):
			start_sequence_position = len(sequence) - len(rc_primer)
			if not nucleotide_comparison(rc_primer[j],sequence[j + start_sequence_position]):
				cpt_bad_r_primers += 1
				result = False
				break

	if result:
		cpt_sequences_matched += 1
		# Writing in the output file	
		# Sequence ID
		output_file.write(fastq_file_lines[i] + "\n")
		# Fasta sequence without primers
		output_file.write(sequence[len(f_primer):start_sequence_position] + "\n")
		# +
		output_file.write('+' + "\n")
		# Quality sequence without primers
		output_file.write(fastq_file_lines[i+3][len(f_primer):start_sequence_position] + "\n")
		seq_length_total += len(sequence[len(f_primer):start_sequence_position])
	else:
		cpt_sequences_unmatched += 1
		# Writing in the output file
		# Sequence ID
		output_file2.write(fastq_file_lines[i] + "\n")
		# Fasta sequence without primers
		output_file2.write(sequence + "\n")
		# +
		output_file2.write('+' + "\n")
		# Quality sequence without primers
		output_file2.write(fastq_file_lines[i+3] + "\n")


output_file.close()
output_file2.close()

if cpt_sequences_matched != 0:
	print "Average length : " + str(seq_length_total/cpt_sequences_matched)
else:
	print "0 Sequence in output file"

print "Number of trimmed sequences : " + str(cpt_sequences_matched)
print "Number of discarded sequences : " + str(cpt_sequences_unmatched)
print "Uncorrect forward primer : " + str(cpt_bad_f_primers) + "\tUncorrect reverse primer : " + str(cpt_bad_r_primers)