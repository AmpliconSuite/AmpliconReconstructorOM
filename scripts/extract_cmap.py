#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

#Extract CMAPs from CMAP file

import argparse


def extract_cmap(cmapF, outname, id_set):
	with open(cmapF) as infile, open(outname,'w') as outfile:
		for line in infile:
			if line.startswith("#"):
				outfile.write(line)
			if not line.startswith("#"):
				fields = line.rsplit("\t")
				if fields[0] in id_set:
					outfile.write(line)


if __name__ == '__main__':
	#Parses the command line arguments
	parser = argparse.ArgumentParser(description="Extract subset of CMAP IDs from CMAP file")
	parser.add_argument("-i", "--input", help="name of input file",required=True)
	parser.add_argument("-l", "--id_list", help="CMAP ID numbers",nargs='+',required=True)
	parser.add_argument("-o", "--outname", help="extracted file name (default to [inputname]_[cmapID1]_..._[cmapIDn].cmap")
	args = parser.parse_args()

	outname = args.outname
	id_set = set(args.id_list)

	if not outname:
		id_string = "_".join(args.id_list)
		outname = args.input[:args.input.index(".cmap")] + "_" + id_string + ".cmap"

	extract_cmap(args.input,outname,id_set)
 
