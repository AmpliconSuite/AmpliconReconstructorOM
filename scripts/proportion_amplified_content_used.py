#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import argparse
import copy
import datetime
import json
import math
import os
from subprocess import call
from Queue import Queue

from intervaltree import IntervalTree

from ContigAlignmentGraph import *
from bionanoUtil import *
from breakpoint_graph import *

# take a graph file and take a cycles file
# report proportion of amplified content explained in each (segment ids). Print missing segment ids (& total length)


def parse_cycles(cyclesf):
    seg_seq_d = {0:(None,0,0)}
    cyclelist = []
    with open(cyclesf) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                seg_seq_d[int(fields[1])] = (fields[2],float(fields[3]),float(fields[4]))

            elif line.startswith("Cycle"):
                fields = line.rstrip().rsplit(";")
                cfd = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                segSeq = [int(x[:-1]) for x in cfd['Segments'].rsplit(',')]
                cyclelist.append(segSeq)

    return seg_seq_d, cyclelist


def graph_amped_segs(graphf, amp_cutoff):
    amped_segs = {}
    unamped_segs = {}
    total_amped = 0.
    total = 0.
    ind = 0
    tot_amp_segs = 0
    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                ind+=1
                fields = line.rsplit()
                chrom = fields[1].rsplit(":")[0]
                p1 = float(fields[1].rsplit(":")[1][:-1])
                p2 = float(fields[2].rsplit(":")[1][:-1])
                cn = float(fields[3])
                size = p2 - p1
                total+=size
                if cn >= amp_cutoff:
                    amped_segs[(chrom, p1, p2)] = ind
                    total_amped+=size
                    tot_amp_segs += 1

                else:
                    unamped_segs[(chrom, p1, p2)] = ind

    print(str(tot_amp_segs) + " of " + str(ind) + " graph segments amplified to at least " + str(amp_cutoff) + " (" +
          str(total_amped) + " bp)")
    print(str(100*total_amped/total) + "% of genomic content in graph amplified")
    return amped_segs, unamped_segs, total_amped


def assess_cycle_amplified_content(amped_segs, unamped_segs, seg_seq_d, cyclelist, total_amped):
    print("CycleID\tN_amp_segs\tN_unamp_segs\tProportionOfGraphAmpExplained\tUnexplainedAmpLength\tUnexplainedSegs")
    for ind, c in enumerate(cyclelist):
        cycle_amped_content = 0.
        cycle_total_content = 0.
        cycle_amp_segs = set()
        cycle_unamp_segs = set()
        for x in set(c):
            if x == 0:
                continue

            ptup = seg_seq_d[x]
            size = ptup[2] - ptup[1]
            cycle_total_content += size
            if ptup in amped_segs:
                cycle_amped_content += size
                cycle_amp_segs.add(ptup)

            elif ptup in unamped_segs:
                cycle_unamp_segs.add(ptup)

            else:
                sys.stderr.write("Graph segment endpoints do not match cycle segment endpoints.")
                sys.exit(1)

        graph_amp_segs = set(amped_segs.keys())
        unused_amp_segs = sorted([str(amped_segs[x]) for x in graph_amp_segs.difference(cycle_amp_segs)])
        unexplained_len = sum([seg_seq_d[int(x)][2] - seg_seq_d[int(x)][1] for x in unused_amp_segs])
        print(str(ind+1) + "\t" + str(len(cycle_amp_segs)) + "\t" + str(len(cycle_unamp_segs)) + "\t" +
              str(cycle_amped_content/total_amped) + "\t" + str(unexplained_len) + "\t" + ",".join(unused_amp_segs))

###"MAIN"###
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="How much amplified content is explained by AR output? "
                                                 "Graph and cycles files must be formatted to identical segments (not"
                                                 "true for direct AA output without running a conversion.)")
    parser.add_argument("-g", help="AA graph file", required=True)
    parser.add_argument("-c", help="AA cycles file", required=True)
    parser.add_argument("--minCN", help="Minimum CN to consider as amplified (default 4.5)", type=float, default=5.0)
    args = parser.parse_args()

    amped_segs, unamped_segs, total_amped = graph_amped_segs(args.g, args.minCN)
    seg_seq_d, cyclelist = parse_cycles(args.c)
    assess_cycle_amplified_content(amped_segs, unamped_segs, seg_seq_d, cyclelist, total_amped)
