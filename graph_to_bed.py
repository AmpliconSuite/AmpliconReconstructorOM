#!/usr/bin/env python

import sys
import os

graph_file = sys.argv[1]
basename = os.path.splitext(sys.argv[1])[0]
with open(graph_file) as infile, open(basename + ".bed", 'w') as outfile:
	for line in infile:
		if line.startswith("sequence"):
			fields = line.rsplit()
			p1,p2 = fields[1],fields[2]
			chrom = p1.rsplit(":")[0]
			start = p1.rsplit(":")[1][:-1]
			end = p2.rsplit(":")[1][:-1]
			outfile.write("\t".join([chrom,start,end]) + "\n")