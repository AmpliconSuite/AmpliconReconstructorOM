#!/usr/bin/env python3

import argparse
import os.path
from collections import defaultdict
import sys

MIN_PYTHON = (3, 0)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

rev_dir = {'+': '-', '-': '+'}


# parse the breakpoint graph to indicate for two endpoints if there is an edge and collect info on the graph segments.
def parse_BPG(BPG_file):
    bidirectional_edge_dict = defaultdict(set)
    seg_end_pos_d = {}
    pos_to_seg = {}
    seqnum = 0
    segToSize = {}
    segToCN = {}
    bps_to_orientations = {}
    with open(BPG_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit()
            if not fields:
                continue

            if fields[0] == "discordant":
                e_rep = fields[1].rsplit("->")
                start = e_rep[0][:-1]
                end = e_rep[1][:-1]
                bidirectional_edge_dict[start].add(end)
                bidirectional_edge_dict[end].add(start)
                bps_to_orientations[(start, end)] = (e_rep[0][-1], e_rep[1][-1])

            elif fields[0] == "sequence":
                seqnum+=1
                seg_end_pos_d[str(seqnum)] = (fields[1][:-1], fields[2][:-1])
                fields = line.rstrip().rsplit()
                lchrom, lpos = fields[1].rsplit(":")
                lpos = int(lpos[:-1])
                rchrom, rpos = fields[2].rsplit(":")
                rpos = int(rpos[:-1])

                cn = float(fields[3])
                segToCN[seqnum] = cn
                segToSize[seqnum] = rpos - lpos
                pos_to_seg[fields[1][:-1]] = (str(seqnum), 'left')
                pos_to_seg[fields[2][:-1]] = (str(seqnum), 'right')

    return bidirectional_edge_dict, seg_end_pos_d, pos_to_seg, segToCN, segToSize, bps_to_orientations


def id_to_loc(BPG_file):
    seqIDD = {}
    with open(BPG_file) as infile:
        seqInd = 0
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                seqInd += 1
                chrom = fields[1].rsplit(":")[0]
                p1 = int(fields[1].rsplit(":")[1][:-1])
                p2 = int(fields[2].rsplit(":")[1][:-1])
                seqIDD[seqInd] = (chrom, p1, p2)

    return seqIDD


def pair_is_adjacent(a_id, b_id, a_dir, b_dir, seqIDD):
    if abs(a_id - b_id) == 1:
        if a_dir == b_dir:
            if seqIDD[a_id][0] != seqIDD[b_id][0]:
                return False

            else:
                if a_dir == "+":
                    if seqIDD[b_id][1] - seqIDD[a_id][2] == 1:
                        return True

                    return False

                else:
                    if seqIDD[a_id][1] - seqIDD[b_id][2] == 1:
                        return True

                    return False

        else:
            return False
    else:
        return False


def pair_is_edge(a_id, b_id, a_dir, b_dir, bpg_dict, seg_end_pos_d):
    try:
        rObj1_end = seg_end_pos_d[a_id][-1] if a_dir == "+" else seg_end_pos_d[a_id][0]
        rObj2_start = seg_end_pos_d[b_id][0] if b_dir == "+" else seg_end_pos_d[b_id][-1]
        return rObj1_end in bpg_dict[rObj2_start], (rObj2_start, rObj1_end)

    except KeyError:
        return False, None


# take a cycles file and a set of cycle IDs
def read_cycles(cfile, elem_to_get):
    with open(cfile) as infile:
        for line in infile:
            if line.startswith("Cycle"):
                fields = line.rstrip().rsplit(";")
                fd = {x.rsplit("=")[0]: x.split("=")[1] for x in fields}
                seq = fd["Segments"]
                cnum = int(fd["Cycle"])
                if cnum == elem_to_get:
                    return seq

    print("Cycle id: " + str(elem_to_get) + ", not found in cycles file")
    return []


def graph_stats(segToSize, segToCN, cutoff=10):
    # full size, amp size, nsegs
    # print(segToCN)
    # print(segToSize)
    nsegs = len(segToSize)
    amp_glen = sum([size for k, size in segToSize.items() if segToCN[k] >= cutoff])
    # print(amp_glen)
    amp_nsegs = len([k for k, cn in segToCN.items() if cn >= cutoff])
    return nsegs, amp_glen, amp_nsegs


def amp_stats(p, sl, amp_glen, amp_nsegs, segToCN, maxInd, cutoff=10):
    num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp = 0, 0, 0, 0
    all_filt_ids = set()
    seq = [int(x[:-1]) for x in p.rsplit(",")]
    # print(seq)
    all_filt_ids |= set([x for x in seq if maxInd >= x > 0])
    num_bp_exp = sum([sl[x] for x in all_filt_ids if segToCN[x] >= cutoff])
    prop_bp_exp = num_bp_exp / amp_glen
    num_gsegs_exp = len(all_filt_ids)
    prop_gsegs_exp = num_gsegs_exp/amp_nsegs

    return num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp


def compare_breakpoints(cp, bpg_dict, seg_end_pos_d, seqIDD, segToCN, cutoff=10):
    tot_amp_AA_bps = 0
    tot_amp_AR_bps = 0
    AA_bps = set()
    amped_AA_bps = set()
    # go through AA breakpoints and check if they are amped
    for bp_s, bp_e_set in bpg_dict.items():
        for bp_e in bp_e_set:
            if tuple(sorted((bp_s, bp_e))) in AA_bps:
                continue

            l = False
            r = False
            for k, et in seg_end_pos_d.items():
                if bp_s == et[0] or bp_s == et[1]:
                    if segToCN[int(k)] >= cutoff:
                        l = True

                if bp_e == et[0] or bp_e == et[1]:
                    if segToCN[int(k)] >= cutoff:
                        r = True

            if l and r:
                amped_AA_bps.add(tuple(sorted((bp_s, bp_e))))
                tot_amp_AA_bps += 1

            AA_bps.add(tuple(sorted((bp_s, bp_e))))
            # AA_bps.add((bp_e, bp_s))

    amplified_AR_bps = set()
    AR_bps = set()
    AR_AA_bps = set()
    oriented_AR_bps = {}
    path = cp.rsplit(",")
    amplified_AR_AA_bps = set()
    if path[0] == "0+":
        path = path[1:-1]
    else:
        path.append(path[0])

    pairs = zip(path[:-1], path[1:])
    for p in pairs:
        # check each pair
        i, i_o = p[0][:-1], p[0][-1]
        t, t_o = p[1][:-1], p[1][-1]
        bpg_adjacency, cbp = pair_is_edge(i, t, i_o, t_o, bpg_dict, seg_end_pos_d)
        ar_e_pair, ar_o_pair = zip(*sorted(zip(cbp, (i_o, t_o))))
        ar_e_pair, ar_o_pair = tuple(ar_e_pair), tuple(ar_o_pair)
        if not pair_is_adjacent(int(i), int(t), i_o, t_o, seqIDD):
            AR_bps.add(ar_e_pair)
            oriented_AR_bps[ar_e_pair] = (ar_o_pair)
            if bpg_adjacency and segToCN[int(i)] > cutoff and segToCN[int(t)] > cutoff:
                AR_AA_bps.add(ar_e_pair)
                amplified_AR_bps.add(ar_e_pair)
                amplified_AR_AA_bps.add(ar_e_pair)

            elif bpg_adjacency:
                # is this connection amped?
                AR_AA_bps.add(ar_e_pair)

            elif segToCN[int(i)] > cutoff and segToCN[int(t)] > cutoff:
                amplified_AR_bps.add(ar_e_pair)
                tot_amp_AR_bps+=1

    # which bps are not used?
    unused_AAs = amped_AA_bps.difference(AR_bps)
    uniq_AR_bps = AR_bps.difference(AR_AA_bps)
    prop_AA_amp_used = len(amplified_AR_AA_bps) / len(amped_AA_bps)
    return tot_amp_AA_bps, amplified_AR_AA_bps, uniq_AR_bps, oriented_AR_bps, prop_AA_amp_used, AA_bps, amped_AA_bps, unused_AAs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Profile an AR reconstruction cycle against AA input data")
    parser.add_argument("-c", "--cycles", help="AA-formatted cycles file CONVERTED TO AA GRAPH COORDINATES",
                        required=True)
    parser.add_argument("--cycle_id", help="Cycle id element to use", type=int, required=True)
    parser.add_argument("-g", "--graph", help="AA-formatted graph file", required=True)
    parser.add_argument("--amp_thresh", help="CN threshold to be considered amplified", type=float, default=10.0)

    args = parser.parse_args()
    elem_to_get = args.cycle_id
    basename = os.path.splitext(os.path.basename(args.cycles))[0]
    oname = basename + "_cycle_" + str(args.cycle_id)

    # retrieve a set of used junctions
    bidirectional_edge_dict, seg_end_pos_d, pos_to_seg, segToCN, segToSize, bps_to_os = parse_BPG(args.graph)
    # retrieve a set of segment IDs (mapped to length and CN)
    nsegs, amp_glen, amp_nsegs = graph_stats(segToSize, segToCN, cutoff=args.amp_thresh)
    # create a map of segment ID to segment coordinates
    seqIDD = id_to_loc(args.graph)
    # retrieve the cycles as an oriented segment list
    cp = read_cycles(args.cycles, elem_to_get)

    # print("Path lengths: ")
    path = cp.rsplit(",")
    if path[0] == "0+":
        path = path[1:-1]

    path_set = set([x[:-1] for x in path])
    sl = sum([segToSize[int(x[:-1])] for x in path])

    # compute statistics about the amplification of segments
    num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp = amp_stats(cp, segToSize, amp_glen, amp_nsegs, segToCN,
                                                                       nsegs, cutoff=args.amp_thresh)

    # compute statistics about the edges used in AR vs AA
    tot_amp_AA_bps, AR_AA_bps, uniq_AR_bps, oriented_AR_bps, bp_amp_frac, AA_bps, amped_AA_bps, unused_AA_bps = \
        compare_breakpoints(cp, bidirectional_edge_dict, seg_end_pos_d, seqIDD, segToCN, cutoff=args.amp_thresh)

    profname = oname + "_profile.tsv"

    # write amplicon profile table
    with open(profname, 'w') as profout:
        profout.write("\nProfile of cycle id:\t" + str(args.cycle_id) + "\n")
        # print("Amplified content used in AA graph vs AR reconstructions")
        profout.write("Amplification threshold selected by user for comparison:\t" + str(args.amp_thresh) + "\n")
        profout.write("Amplified AA graph genomic content (bp):\t" + str(amp_glen) + "\n")
        profout.write("Amplified AA graph segments:\t" + str(amp_nsegs) + "\n")
        profout.write("Amplified AA graph breakpoints:\t" + str(tot_amp_AA_bps) + "\n")

        profout.write("Length of reconstructed cycleq (bp):\t" + str(sl) + "\n")
        profout.write("Number of graph segments used in cycle:\t" + str(len(path_set)) + "\n")

        profout.write("Amplified genomic content from graph used in cycle (bp)|(%):\t" + str(num_bp_exp) + "\t" +
                      '{:.1%}'.format(prop_bp_exp) + "\n")

        profout.write("Amplified segments from graph used in cycle (number)|(%):\t" + str(num_gsegs_exp) + "\t"
                      '{:.1%}'.format(prop_gsegs_exp) + "\n")

        # TODO
        profout.write("Amplified graph breakpoints used in cycle (number)|(%):\t" + str(len(AR_AA_bps)) + "\t" +
                      '{:.1%}'.format(bp_amp_frac) + "\n")
        profout.write("Number of reconstruction breakpoints not in graph:\t" + str(len(uniq_AR_bps)) + "\n")

        # profout.write("Median estimated CN of amplicon (controlled for duplicated segments):\t" + "XXX" + "\n")

    # write segment usage table
    seg2PathMultiplicity = defaultdict(int)
    for x in path:
        cid = int(x[:-1])
        seg2PathMultiplicity[cid]+=1

    with open(oname + "_segments_profile.tsv", 'w') as outfile:
        outfile.write("chrom\tstart\tend\tsegment_id\tCN\tmultiplicity_in_reconstruction\n")
        for cid in range(1, nsegs+1):
            a, b = seg_end_pos_d[str(cid)]
            chrom, p1 = a.rsplit(":")
            _, p2 = b.rsplit(":")
            cn = segToCN[cid]
            outfile.write("\t".join([str(cid), chrom, p1, p2, str(cid), str(cn), str(seg2PathMultiplicity[cid])]) + "\n")

    # write breakpoint usage table
    with open(oname + "_breakpoint_profile.tsv", 'w') as outfile:
        outfile.write("chrom1\tpos1\tchrom2\tpos2\torientation1\torientation2\tamplified\tin_reconstruction\tjunction_not_in_graph\n")
        for x in AA_bps:
            if x in bps_to_os:
                bpos = bps_to_os[x]
            else:
                bpos = bps_to_os[(x[1], x[0])]
                bpos = bpos[::-1]

            bpos = [bpos[0], rev_dir[bpos[1]]]

            isAmp = "True" if x in amped_AA_bps else "False"
            inRecons = "True" if x in AR_AA_bps else "False"
            olist = x[0].rsplit(":") + x[1].rsplit(":") + bpos + [isAmp, inRecons, "False"]
            outfile.write("\t".join(olist) + "\n")

        for y in uniq_AR_bps:
            bpos = list(oriented_AR_bps[y])
            s1 = segToCN[int(pos_to_seg[y[0]][0])]
            s2 = segToCN[int(pos_to_seg[y[1]][0])]
            isAmp = s1 >= args.amp_thresh and s2 >= args.amp_thresh
            olist = y[0].rsplit(":") + y[1].rsplit(":") + bpos + [str(isAmp), "True", "True"]
            outfile.write("\t".join(olist) + "\n")
