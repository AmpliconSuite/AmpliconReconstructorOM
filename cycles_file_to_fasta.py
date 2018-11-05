#!/usr/bin/python
import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def readRefG(GENOME_BUILD_FASTA):
	#read human ref genome
	rgLength = 0
	chrLengthD = {}
	seqD = {}
	sys.stdout.write("READING REF GENOME FASTA\n")
	fasta_sequences = SeqIO.parse(open(GENOME_BUILD_FASTA),'fasta')
	for fasta in fasta_sequences:
		name = fasta.id
		#if name == "HPV16REF.1":
			#name = "hpv16ref_1"
		#sys.stdout.write("FINISHED READING " + name + "\n")
		sequence = str(fasta.seq)
		seqD[name] = sequence
		rgLength+=len(sequence)
		print name, len(sequence)

		chrLengthD[name] = len(sequence)

	return rgLength, chrLengthD, seqD

def reconstructEpisomeCycle(bpointf,outpre):
	#rev comp
	#str(Seq(safeSeg).reverse_complement())
	cycleFasta = open(outpre + "_cycles.fasta",'w')
	infile = open(bpointf)
	for line in infile:
		if "Cycle=" in line:
			fields = line.rstrip().rsplit(";")
			lineD = {x.rsplit("=")[0]:x.rsplit("=")[1] for x in fields}
			segs = lineD["Segments"].rsplit(",")
			recSeq = ""
			for i in segs:
				seg = i[:-1]
				if seg != "0":
					strand = i[-1]
					segSeq = segSeqD[seg]
					if strand == "-":
						segSeq = str(Seq(segSeq).reverse_complement())
					recSeq+=segSeq

			outname = args.input.rsplit(".")[0]
			outname = outname + "_cycle_" + lineD["Cycle"]
			cycleFasta.write(">" + outname + "\n")
			cycleFasta.write(recSeq + "\n")

	cycleFasta.close()

def segsToSeq(bpointf,outpre):
	segFasta = open(outpre + "_segments.fasta",'w')
	segSeqD = {}
	infile = open(bpointf)
	for line in infile:
		if line.startswith("Segment"):
			fields = line.rstrip().rsplit()
			lowerBound = int(fields[3])
			upperBound = int(fields[4])
			relSeq = seqD[fields[2]][lowerBound:upperBound]
			segFasta.write(">" + "_".join(fields) + "\n")
			segFasta.write(relSeq + "\n")
			segNum = fields[1]
			segSeqD[segNum] = relSeq

	segFasta.close()
	return segSeqD


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="cyle text file",required=True)
parser.add_argument("-r", "--ref", help="reference genome",required=True)
parser.add_argument("-o", "--output", help="output prefix",required=True)
args = parser.parse_args()

rgLength,chrLength,seqD = readRefG(args.ref)
segSeqD = segsToSeq(args.input,args.output)
reconstructEpisomeCycle(args.input,args.output)
