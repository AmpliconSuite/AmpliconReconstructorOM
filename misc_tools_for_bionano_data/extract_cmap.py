#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#Extract CMAP from CMAP file

import sys
import os
import argparse

def extract_cmap(cmapF,outname,id):
	with open(cmapF) as infile, open(outname,'w') as outfile:
		for line in infile:
			if line.startswith("#"):
				outfile.write(line)
			if not line.startswith("#"):
				fields = line.rsplit("\t")
				if fields[0] == id:
					outfile.write(line)


if __name__ == '__main__':
	#Parses the command line arguments
	parser = argparse.ArgumentParser(description="Extract CMAP from CMAP file")
	parser.add_argument("-i", "--input", help="name of input file",required=True)
	parser.add_argument("-m", "--mapID", help="CMAP ID number",required=True)	
	parser.add_argument("-o", "--outname", help="extracted file name (default to [inputname]_[cmapID].cmap")
	args = parser.parse_args()

	outname = args.outname
	if not outname:
		#infName = args.input.rsplit("/")[-1].rsplit(".")[0]
		#outname = "/".join(args.input.rsplit("/")[:-1]) + "/" + infName + "_CMAPID" + args.mapID + ".cmap"
		outname = args.input[:args.input.index(".cmap")] + "_CMAPID" + args.mapID + args.input[args.input.index(".cmap"):]

	extract_cmap(args.input,outname,args.mapID)
 