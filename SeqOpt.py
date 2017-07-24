#!/usr/bin/env python

## ----------------------------------------
## Sequence Optimzer
## ----------------------------------------
## Quickly modify codon usage in a given
## sequence subject to provided weights
## ----------------------------------------
## Author:  Aditya Radhakrishnan
## Contact: adityaradhakrishnan@gmail.com
## Website: radhakrishnan.me
## Date:    2017/07/23
## ----------------------------------------

import sys

def FASTAParser(FileName = ''):

	try:
		FASTAFile = open(FileName, 'r')
	except IOError:
		sys.exit('FASTA file not found at specified location. Check name and directory!')

	# Extract name of sequence from first line of FASTAFile, accounting for '>' and '\n'

	SeqName = FASTAFile.readline()[1:-1]

	# Store sequence of gene as list of codons. Use rstrip instead of [:-1] in case last line has
	# potential issues with inclusion of the newline character.

	Sequence = ""

	for iLine in FASTAFile:
		Sequence += iLine.rstrip("\n")

	SeqSplit = [Sequence[iX:iX+3] for iX in range(0, len(Sequence), 3)]

	# Check that the sequence is actually proper

	if len(SeqSplit[-1]) != 3:
		sys.exit('This sequence doesn\'t encode a protein. It\'s not a multiple of 3 in length!')
	elif SeqSplit[0] != 'ATG':
		sys.exit('This sequence doesnt\'t start with an ATG! Check your FASTA File!')

	# Return name and split sequence

	return [SeqName, SeqSplit]


[FASTAName, FASTASplit] = FASTAParser("Sample.fasta")
