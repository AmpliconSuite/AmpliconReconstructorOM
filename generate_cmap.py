#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import sys
import os
import re
import argparse
from itertools import groupby
from collections import defaultdict

complementary_nucleotide = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


"""
Modified method from brentp on BioStars
given a fasta file. yield dictionary of header, sequence
"""
def fasta_reader(fasta_file, chroms_to_get, getAll=False):
    fasta_dict = {}
    with open(fasta_file) as infile:
        faiter = (x[1] for x in groupby(infile, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            seq_name = next(header)[1:].rstrip().rsplit()[0]
            if (seq_name in chroms_to_get) or getAll:
                print("Reading " + seq_name)
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in next(faiter))
                fasta_dict[seq_name] = seq

    return fasta_dict


def rev_complement(seq):
    return ''.join([complementary_nucleotide[a] for a in seq[::-1]])


def read_graph(graphF):
    with open(graphF) as infile:
        segSeqL = []
        chroms_to_get = set()
        head = next(infile).rstrip().split()
        segN = 0
        for line in infile:
            if line.startswith("sequence"):
                segN += 1
                fields = line.rstrip().split()
                f1 = fields[1].rsplit(":")
                f2 = fields[2].rsplit(":")
                chrom = f1[0]
                chroms_to_get.add(chrom)
                p1 = int(f1[1][:-1])
                p2 = int(f2[1][:-1])
                lowerBound = min(p1, p2)
                upperBound = max(p1, p2)
                segSeqL.append((fields[1] + "|" + fields[2], chrom, lowerBound - 1, upperBound - 1))

    return segSeqL, chroms_to_get


def read_bed(bedfile):
    bed_dict = defaultdict(list)
    with open(bedfile) as infile:
        for line in infile:
            if not line.startswith("#"):
                fields = line.rstrip().rsplit()
                bed_dict[fields[0]].append((int(fields[1]), int(fields[2])))

    chroms_to_get = set(bed_dict.keys())
    segSeqL = []
    for chrom, region_list in bed_dict.items():
        for p in region_list:
            pstring = str(p[0]) + "-|" + str(p[1]) + "+"
            segSeqL.append((pstring, chrom, p[0] - 1, p[1] - 1))

    return segSeqL, chroms_to_get


def segsToSeq(segSeqL, seqD):
    segSeqD = {}
    for i in segSeqL:
        segSeqD[i[0]] = seqD[i[1]][i[2]:i[3]]

    return segSeqD


def makeCMAP(outprefix, segSeqL, segSeqD, enzyme, regExpEnzTup, minLabel, minSize):
    print("Generating CMAP for enzyme " + enzyme)
    cmap_header = "# CMAP File Version:	0.1\n# Label Channels:	1\n# Nickase Recognition Site 1:	" + enzyme + "\n# Enzyme1:	" + enzyme + "\n# Number of Consensus Nanomaps:	" + str(
        len(segSeqD)) + "\n#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence\n#f int	float	int	int	int	float	float	int	int\n"
    key_header = "# CMAP = " + args.output + "_key.txt\n# filter: Minimum Labels = " + str(
        minLabel) + "\n# filter: Minimum Size (Kb) = " + str(minSize / 1000.0) + "\n"
    regExpEnz, offset = regExpEnzTup
    # cmap 1 based
    offset += 1
    with open(outprefix + ".cmap", 'w') as cmap_outfile, open(outprefix + "_key.txt", 'w') as key_outfile:
        cmap_outfile.write(cmap_header)
        key_outfile.write(key_header)
        for CMapNum, i_tup in enumerate(segSeqL):
            i = i_tup[0]
            CMapId = str(CMapNum + 1)
            currSeq = segSeqD[i].upper()
            mapLen = str(float(len(currSeq)))
            if mapLen < minSize:
                continue

            print("Generating map for " + i)
            revRegExpEnz = rev_complement(regExpEnz)
            rev_offset = len(regExpEnz) - offset
            indVPos = [float(m.start() + offset) for m in re.finditer('(?=' + regExpEnz + ')', currSeq)]
            indVNeg = [float(m.start()) + rev_offset for m in re.finditer('(?=' + revRegExpEnz + ')', currSeq)]
            # allPosns = merge(indVPos,indVNeg) #for 1,000,000 total sites or more this would be better
            allPosns = sorted(list(set(indVNeg + indVPos)))
            if len(allPosns) < minLabel:
                continue

            totalLabels = str(len(allPosns))
            # CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
            for ind, val in enumerate(allPosns):
                cmap_outfile.write("\t".join([CMapId, mapLen, totalLabels, str(ind + 1), "1", str(val), "1.0\t1\t1\n"]))

            cmap_outfile.write(
                "\t".join([CMapId, mapLen, totalLabels, str(len(allPosns) + 1), "0", mapLen, "1.0\t1\t1\n"]))
            key_outfile.write("\t".join([CMapId, i, mapLen]) + "\n")

    print("Finished")


parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", help="reference genome fasta", required=True)
parser.add_argument("-o", "--output", help="output prefix")
parser.add_argument("-e", "--enzyme", help="restriction enzyme: BspQI, BbvCI, BsmI, BsrDI, DLE1", required=True)
parser.add_argument("-l", "--labels", type=int, help="minimum number of labels in reported CMAP, default 0", default=0)
parser.add_argument("-s", "--size", type=float, help="minimum map size for reported CMAP (basepairs), default 0",
                    default=0)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-g", "--graph",
                   help="breakpoint graph .txt file, if not supplied, --makeRef or --bed must be supplied")
group.add_argument("--makeRef", help="Make CMAP from entire reference fasta", action='store_true')
group.add_argument("-b", "--bed",
                   help="breakpoint graph bed file, it not supplied, --makeRef or --graph must be supplied")

args = parser.parse_args()

enzymeSequences = {
    "BSPQI": ("GCTCTTC", 8),
    "BBVCI": ("CCTCAGC", 2),
    "BSMI": ("GAATGC", 4),
    "BSRDI": ("GCAATG", 5),
    "DLE1": ("CTTAAG", 2)
}

# "bseCI":("ATCGAT",),
# "DCM1":"CCAGG",
# "DCM2":"CCTGG"

if args.labels < 0:
    sys.exit("Innaproriate min number of labels " + str(args.labels) + ", exiting")

if args.size < 0:
    sys.exit("Innaproriate min map size " + str(args.size) + ", exiting")

args.enzyme = args.enzyme.upper()
if args.enzyme not in enzymeSequences:
    print("Valid enzymes are " + str(enzymeSequences.keys()))
    sys.exit("Unrecognized enzyme " + args.enzyme + ", exiting")

if not args.output:
    if args.graph:
        args.output = ".".join(args.graph.rsplit(".")[:-1])

    elif args.bed:
        args.output = ".".join(args.bed.rsplit(".")[:-1])

    else:
        args.output = ".".join(args.ref.rsplit("/")[-1].rsplit(".")[:-1])

outprefix = args.output
if args.enzyme not in args.output:
    outprefix = outprefix + "_" + args.enzyme

regExpEnzTup = enzymeSequences[args.enzyme]

if args.graph:
    segSeqL, chroms_to_get = read_graph(args.graph)
    seqD = fasta_reader(args.ref, chroms_to_get)

elif args.bed:
    segSeqL, chroms_to_get = read_bed(args.bed)
    seqD = fasta_reader(args.ref, chroms_to_get)

else:
    chroms_to_get = set()
    seqD = fasta_reader(args.ref, chroms_to_get, True)
    segSeqL = []
    for i in sorted(seqD, key=lambda x: x.lstrip('chr')):
        seqLen = len(seqD[i])
        if seqLen > args.size:
            segSeqL.append((i, i, 0, len(seqD[i])))

segSeqD = segsToSeq(segSeqL, seqD)
makeCMAP(outprefix, segSeqL, segSeqD, args.enzyme, regExpEnzTup, args.labels, args.size)
