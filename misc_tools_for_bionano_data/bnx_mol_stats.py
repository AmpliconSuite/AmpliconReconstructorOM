#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import sys
import os
import argparse
from matplotlib import pyplot as plt
import numpy as np


#parses bnx file and reports molecule stats (rudimentary)


def parseBNXFile(fName): 
	moleculeD = {}
	with open(fName) as infile:
		for line in infile:
			if line.startswith("#0h"):
					head = line.rstrip().rsplit()[1:]
			elif line.startswith("0"):
				fields = line.rstrip().rsplit()
				#print head
				#print fields
				fD = dict(zip(head,fields))
				#get the chrom
				moleculeD[fD["MoleculeID"]] = float(fD["Length"])

	return moleculeD

#plot it out
def plot_lens(molD):
	lens = molD.values()
	print np.mean(lens),np.std(lens)
	fig = plt.figure()
	plt.hist(lens,bins=40)
	plt.savefig(outName + "_mol_len_hist.png",dpi=300)
	plt.close()



#inputs must be 
#-contig name
#-map file
#-all_bnx file

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract aligned bionano molecules from contigs \
    	at given reference genome positions")
	# parser.add_argument("-m","--map",help="pileup map file",required=True)
	# parser.add_argument("-c","--contig",help="contig ID to get molecules from",required=True)
	# parser.add_argument("--coords",help="input chromosome cordinates, formatted as 'start-end', \
 #    	default entire chromosome")
	parser.add_argument("--bnx",help="merged BNX file",required=True)
	parser.add_argument("-o","--outname",help="output filename prefix",default="")
	args = parser.parse_args()

	if args.outname:
		args.outname+="_"
	outName = args.outname

	molD = parseBNXFile(args.bnx)

	plot_lens(molD)