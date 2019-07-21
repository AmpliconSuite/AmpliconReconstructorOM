#!/usr/bin/env python

import sys
import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Corrects and extends alignment paths to produce BioNano contig/AA segment scaffolds")
parser.add_argument("-g","--graph",help="AA-formatted graph file",required=True)
parser.add_argument("--padding",type=int,
	help="padding to add on either side of graph regions, in bp. Used in creation of padded_sorted file",required=True, default=20000)
args = parser.parse_args()

graph_file = args.graph
padding = args.padding
basename = os.path.splitext(args.graph)[0]
bed_dict = defaultdict(list)
with open(graph_file) as infile, open(basename + ".bed", 'w') as outfile:
	for line in infile:
		if line.startswith("sequence"):
			fields = line.rsplit()
			p1,p2 = fields[1],fields[2]
			chrom = p1.rsplit(":")[0]
			start = p1.rsplit(":")[1][:-1]
			end = p2.rsplit(":")[1][:-1]
			bed_dict[chrom].append((int(start),int(end)))
			outfile.write("\t".join([chrom,start,end]) + "\n")

with open(basename + "_padded_sorted_merged.bed",'w') as outfile:
	for chrom in sorted(bed_dict.keys()):
		padded_pos_list = [(max(0,x[0]-padding),x[1]+padding) for x in bed_dict[chrom]]
		sorted_interval_list = sorted(padded_pos_list)
		curr_start,curr_end = sorted_interval_list[0]
		ind = 1
		while ind != len(sorted_interval_list):
			s,e = sorted_interval_list[ind]
			if s > curr_end+1:
				outfile.write("\t".join([chrom,str(curr_start),str(curr_end)]) + "\n")
				curr_start,curr_end = s,e

			else:
				curr_end = max(e,curr_end)

			ind+=1

		outfile.write("\t".join([chrom,str(curr_start),str(curr_end)]) + "\n")



