#!/usr/bin/env python

import os
import sys

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

cmapfile = sys.argv[1] # seg cmap file
alignment_file = sys.argv[2] # modified to contain the alignment path

rdd = {'+':'-', '-':'+'}

curr_end = 0
idlens = {}
with open(cmapfile) as infile:
    for line in infile:
        if line.startswith("#"):
            continue

        fields = line.rstrip().rsplit()
        cmid = fields[0]
        idlens[cmid] = int(fields[2])

oname = os.path.splitext(os.path.basename(alignment_file))[0]
oname += "_remapped_seg_aln.txt"

with open(alignment_file) as infile, open(oname, 'w') as ofile:
    flip = "_r_" in alignment_file
    ofile.write(next(infile))
    seqline = next(infile)
    ofile.write(seqline)
    ahead = next(infile)
    ahead_fields = ahead.rsplit("\t")
    ahead_fields = ahead_fields[:7] + ["imputed",] + ahead_fields[7:]
    ahead_new = "\t".join(ahead_fields)
    ofile.write(ahead_new)
    seqline = seqline[1:].rsplit("\t")[0]
    seq = [(x[:-1], x[-1]) for x in seqline.rsplit(",")]
    if seq[0][0] == '0':
        seq = seq[1:-1]
    l2seq = {}
    l2sl = {}
    l2an = {}
    l2dir = {}
    an = -1
    finish = 0
    cumstart = 1
    for s, dir in seq:
        an += 1
        nl = idlens[s]
        if nl > 0:
            finish = cumstart+nl

        # print(cumstart,finish)
        for ind, i in enumerate(range(cumstart, finish)):
            l2seq[i] = s
            if dir == "+":
                l2sl[i] = ind+1
            else:
                l2sl[i] = nl - ind
            l2an[i] = an
            l2dir[i] = dir

        if nl > 0:
            cumstart = finish

    # print(l2seq)

    for line in infile:
        fields = line.rstrip().rsplit("\t")
        cs = str(l2seq[int(fields[3])])
        cl = str(l2sl[int(fields[3])])
        segdir = str(l2dir[int(fields[3])])
        if flip:
            segdir = rdd[segdir]

        can = str(l2an[int(fields[3])])
        oline = "\t".join([fields[0], cs, fields[2], cl, '+', segdir, can, "0\t0\t0"]) + "\n"
        ofile.write(oline)
