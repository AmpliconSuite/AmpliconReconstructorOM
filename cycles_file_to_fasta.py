#!/usr/bin/env python
import sys
import os
import argparse

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

"""
Modified method from brentp on BioStars
given a fasta file. yield dictionary of header, sequence
"""
def fasta_reader(fasta_file,chroms_to_get,getAll = False):
	fasta_dict = {}
	with open(fasta_file) as infile:
		faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
		for header in faiter:
			# drop the ">"
			seq_name = header.next()[1:].strip()
			if (seq_name in chroms_to_get) or getAll:
				print("Reading " + seq_name)
				# join all sequence lines to one.
				seq = "".join(s.strip() for s in faiter.next())
				fasta_dict[seq_name] = seq

	return fasta_dict

def writeCycleFasta(cyclef,outpre):
	with open(outpre + "_cycles.fasta",'w') as cycleFasta, open(cyclef) as infile:
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

#get the relevant chromosomes the cycles file
def relChroms(cyclef):
	chromSet = set()
	with open(cyclef) as infile:
		for line in infile:
			if line.startswith("Segment"):
				fields = line.rstrip().rsplit()
				chromSet.add(fields[2])

	return chromSet


def segsToSeq(cyclef,outpre):
	segSeqD = {}
	with open(cyclef) as infile, open(outpre + "_segments.fasta",'w') as segFasta:
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

	return segSeqD


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cycles", help="AA-formatted cycles file",required=True)
parser.add_argument("-r", "--ref", help="path to reference genome",required=True)
parser.add_argument("-o", "--output", help="output prefix",required=True)
args = parser.parse_args()

relChromSet = relChroms(args.cycles)
seqD = fasta_reader(args.ref,relChromSet)
segSeqD = segsToSeq(args.cycles,args.output)
writeCycleFasta(args.cycles,args.output)
