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
import math

from random import shuffle

def FASTAErrorHandler(Split, SName, FName):

	if len(Split[-1]) != 3:
		sys.exit('Error in ' + FName + ' in the sequence ' + SName + '. This sequence doesn\'t encode a protein. It\'s not a multiple of 3 in length!')
	elif Split[0] != 'ATG':
		sys.exit('Error in ' + FName + ' in the sequence ' + SName + '. This sequence doesnt\'t start with an ATG! Check your FASTA File!')



def FASTAParser(FileName = ''):

	try:
		FASTAFile = open(FileName, 'r')
	except IOError:
		sys.exit('FASTA file not found at specified location. Check name & directory!')

	# Extract name of sequence from first line of FASTAFile, accounting for '>' and '\n'

	SeqName = FASTAFile.readline()[1:-1]

	# Store sequence of gene as list of codons. Use rstrip instead of [:-1] in case last line has
	# potential issues with inclusion of the newline character.

	Sequence = ""

	for iLine in FASTAFile:
		Sequence += iLine.rstrip("\n")

	SeqSplit = [Sequence[iX:iX+3] for iX in range(0, len(Sequence), 3)]

	# Check that the sequence is actually proper

	FASTAErrorHandler(SeqSplit, SeqName, FileName)

	# Return name and split sequence

	return [SeqName, SeqSplit]

def OptimalityParser(FileName = '', SynCodonDict = {}):

	# Get codon specific weights for codon adaptation index. Supply a list of well expressed
	# ORFs in your organism of choice to extract the relevant values.

	try:
		OptFile = open(FileName, 'r')
	except IOError:
		sys.exit('Optimized FASTA file not found at specified location. Check name & directory!')

	# We don't really care about the names of these genes, so just extract the sequence info:

	Sequence = ""
	SeqName  = ""
	OptDict  = {}

	for iLine in OptFile:
		if iLine[0] == ">":
			pass
		else:
			Sequence += iLine.rstrip("\n")

	SeqSplit = [Sequence[iX:iX+3] for iX in range(0, len(Sequence), 3)]
	FASTAErrorHandler(SeqSplit, SeqName, FileName)

	for iCodon in SeqSplit:
		if iCodon in OptDict:
			OptDict[iCodon] += 1
		else:
			OptDict[iCodon] = 0.0


	for iAA in SynCodonDict:
		maxAA    = 0
		maxCodon = ''

		for iCodon in SynCodonDict[iAA]:
			try:
				OptDict[iCodon]
			except KeyError:
				sys.exit('You need more sequences in your reference set! The codon ' + iCodon + ' didn\'t even appear once in the reference set')

			if OptDict[iCodon] > maxAA:
				maxAA    = OptDict[iCodon]
				maxCodon = iCodon

		for iCodon in SynCodonDict[iAA]:
			OptDict[iCodon] = OptDict[iCodon]/maxAA

	return OptDict

def SeqSlicer(ProteinSeq):
	print ""
	Range = raw_input("Enter the AA range you want to (de)-optimize (e.g. 96-137): ")
	Index = [int(Idx) - 1 for Idx in Range.split("-")]

	print ""
	print "To verify, the sequence you want to modify is:"

	print ProteinSeq[Index[0]:Index[1]]

	print ""
	Response = raw_input("Is this correct? (Y/N): ")

	if Response == "Y":
		return [Index[0], Index[1], ProteinSeq[Index[0]:Index[1]]]
	else:
		return [Index[0], Index[1], ""]

def OptimizeSequence(SplitSeq, SynCodonDict, CodonTabDict, OptDict):

	ProtSeq = ""
	SubSeq  = [0, 0, ""]

	for iCodon in SplitSeq:
		ProtSeq += CodonTabDict[iCodon]

	print "Your DNA sequence encodes the protein:"
	print ProtSeq
	
	while SubSeq[2] == "":
		SubSeq = SeqSlicer(ProtSeq)

	Length = SubSeq[1] - SubSeq[0]
	CAI    = 0.0

	for iX in range(SubSeq[0], SubSeq[1]):
		CAI += math.log10(OptDict[SplitSeq[iX]]) 

	CAI = 10**(CAI/(Length))

	print ""
	print "The optimality of this subregion as quantified by CAI is: "
	print CAI

	print ""
	Target  = input("What is the target CAI for this region?: ")
	Target  = float(Target)
	NewSeq  = [""]*(Length)
	Indices = range(SubSeq[0], SubSeq[1])

	# Randomly perturb codons so you don't lead to waves or neighbor affects

	shuffle(Indices)

	for iX in Indices:
		AA       = CodonTabDict[SplitSeq[iX]]
		Curr     = SplitSeq[iX]

		for iCodon in SynCodonDict[AA]:
			Delta = math.log10(CAI)*Length - math.log10(OptDict[Curr]) + math.log10(OptDict[iCodon])
			Delta = 10**(Delta/Length)

			if math.fabs(Target - Delta) < math.fabs(Target - CAI):
				Curr = iCodon
				CAI  = Delta

		NewSeq[iX - SubSeq[0]] = Curr

	print "The codon modified sequence is now: "
	print "".join(NewSeq)

	CAI = 0.0

	for iX in NewSeq:
		CAI += math.log10(OptDict[iX]) 

	CAI = 10**(CAI/len(NewSeq))

	print ""
	print "The optimality of this sequence is now:"
	print CAI

SynonymousCodons = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'],
                    'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                    'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
                    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                    'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], '*': ['TAA', 'TAG', 'TGA']
                    }

CodonTable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        	  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        	  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        	  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        	  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        	  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        	  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        	  'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        	  }

[FASTAName, FASTASplit] = FASTAParser("FASTA/Sample.fasta")
OptimalityDict          = OptimalityParser("FASTA/Optimized.fasta", SynonymousCodons)

for iX in OptimalityDict:
	print iX, OptimalityDict[iX]

OptimizeSequence(FASTASplit, SynonymousCodons, CodonTable, OptimalityDict)


