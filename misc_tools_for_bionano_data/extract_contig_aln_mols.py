#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import sys
import os
import argparse

#parses map file and extracts molecule ids aligning to particular coordinates
def parseMapFile(fName,contID,start,end): 
	moleculeS = set()
	head = []
	with open(fName) as infile:
		for line in infile:
			if not line.startswith("#"):
				if not head:
					head = line.rstrip().rsplit("\t")
				elif line:
					fields = line.rstrip().rsplit("\t")
					reducedHead = head[:len(fields)]
					lineD = dict(zip(reducedHead,fields))
					#get the chrom
					if lineD["ContigId"] == contID: 
						#get the posns
						beginM,endM = int(lineD["StartMatchLocation"]),int(lineD["EndMatchLocation"])
						#compare window
						if (beginM < start and endM > start) or (beginM > start and beginM < end):
							moleculeS.add(lineD["MoleculeId"])

	return moleculeS

#extacts molecules from an existing file and prints them out
def writeBNX(bnx_fname,moleculeS,outName):
	with open(bnx_fname) as infile, open(outName,'w') as outfile:
		for line in infile:
			if line.startswith("#"):
				outfile.write(line)
			else:
				if line.startswith("0"):
					printEntry = False
					fields = line.rstrip().rsplit("\t")
					molID = fields[1]
					if molID in moleculeS:
						printEntry = True

				if printEntry:
					outfile.write(line)


#inputs must be 
#-contig name
#-map file
#-all_bnx file

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract aligned bionano molecules from contigs \
    	at given reference genome positions")
	parser.add_argument("-m","--map",help="pileup map file",required=True)
	parser.add_argument("-c","--contig",help="contig ID to get molecules from",required=True)
	parser.add_argument("--coords",help="input chromosome cordinates, formatted as 'start-end', \
    	default entire chromosome")
	parser.add_argument("--bnx",help="merged BNX file",required=True)
	parser.add_argument("-o","--outname",help="output filename prefix",default="")
	args = parser.parse_args()

	if not args.coords:
		start = 0
		end = sys.maxint

	else:
		start = int(args.coords.rsplit("-")[0])
		end = int(args.coords.rsplit("-")[1])

	print "searching for coords " + str(start) + " " + str(end) 
	if end == sys.maxint:
		region = "ALL"
	else:
		region = args.coords

	if args.outname:
		args.outname+="_"
	outName = args.outname + "extracted_molecules_" + args.contig + "_" + region + ".bnx"

	print "Parsing .MAP"
	moleculeS = parseMapFile(args.map,args.contig,start,end)
	print str(len(moleculeS)) + " molecules identified aligning to query region"
	print "Writing BNX"
	writeBNX(args.bnx,moleculeS,outName)