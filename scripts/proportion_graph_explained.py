#!/usr/bin/env python3

import argparse
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
    seqnum = 0
    segToSize = {}
    segToCN = {}
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

    return bidirectional_edge_dict, seg_end_pos_d, segToCN, segToSize


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


def read_cycles(cfile, elems_to_get):
    seqlist = []
    with open(cfile) as infile:
        for line in infile:
            if line.startswith("Cycle"):
                fields = line.rstrip().rsplit(";")
                fd = {x.rsplit("=")[0]: x.split("=")[1] for x in fields}
                seq = fd["Segments"]
                cnum = int(fd["Cycle"])
                if cnum in elems_to_get or not elems_to_get:
                    seqlist.append(seq)

    return seqlist


def graph_stats(segToSize, segToCN, cutoff=10):
    # full size, amp size, nsegs
    nsegs = len(segToSize)
    amp_glen = sum([size for k, size in segToSize.items() if segToCN[k] >= cutoff])
    amp_nsegs = len([k for k, cn in segToCN.items() if cn >= cutoff])
    return nsegs, amp_glen, amp_nsegs


def amp_stats(recons_paths, sl, amp_glen, amp_nsegs, segToCN, maxInd, cutoff=10):
    num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp = 0, 0, 0, 0
    all_filt_ids = set()
    for ind, p in enumerate(recons_paths):
        seq = [int(x[:-1]) for x in p.rsplit(",")]
        all_filt_ids |= set([x for x in seq if maxInd >= x > 0])

    num_bp_exp = sum([sl[x] for x in all_filt_ids if segToCN[x] >= cutoff])
    prop_bp_exp = num_bp_exp / amp_glen
    num_gsegs_exp = len(all_filt_ids)
    prop_gsegs_exp = num_gsegs_exp/amp_nsegs

    return num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp


def compare_breakpoints(recons_paths, bpg_dict, seg_end_pos_d, seqIDD, segToCN, cutoff=10):
    tot_amp_AA_bps = 0
    tot_amp_AR_bps = 0
    used_AA_bps = set()
    used_amped_AA_bps = set()
    # go through AA breakpoints and check if they are amped
    for bp_s, bp_e_set in bpg_dict.items():
        for bp_e in bp_e_set:
            if (bp_s, bp_e) in used_AA_bps:
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
                used_amped_AA_bps.add((bp_s, bp_e))
                tot_amp_AA_bps += 1

            used_AA_bps.add((bp_s, bp_e))
            used_AA_bps.add((bp_e, bp_s))

    used_AR_bps = set()
    used_AR_bps_AA_form = set()
    for cp in recons_paths:
        path = cp.rsplit(",")
        if path[0] == "0+":
            path = path[1:-1]

        if len(path) < 2:
            continue

        pairs = zip(path[:-1], path[1:])
        for p in pairs:
            # check each pair
            i, i_o = p[0][:-1], p[0][-1]
            t, t_o = p[1][:-1], p[1][-1]
            bpg_adjacency, cbp = pair_is_edge(i, t, i_o, t_o, bpg_dict, seg_end_pos_d)
            used_AR_bps_AA_form.add(cbp)
            used_AR_bps_AA_form.add((cbp[1], cbp[0]))
            if not pair_is_adjacent(int(i), int(t), i_o, t_o, seqIDD):
                if bpg_adjacency:
                    fp = ((i, i_o), (t, t_o))
                    rp = ((t, rev_dir[t_o]), (i, rev_dir[i_o]))
                    if fp not in used_AR_bps and rp not in used_AR_bps:
                        used_AR_bps.add(fp)
                        used_AR_bps.add(rp)
                        # is this connection amped?
                        if segToCN[int(i)] > cutoff and segToCN[int(t)] > cutoff:
                            tot_amp_AR_bps+=1

    # which bps are not used?
    unused_AAs = used_amped_AA_bps - used_AR_bps_AA_form
    return tot_amp_AA_bps, tot_amp_AR_bps, tot_amp_AR_bps/tot_amp_AA_bps, unused_AAs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify AA amplicon type")
    parser.add_argument("-c", "--cycles", help="AA-formatted cycles file CONVERTED TO AA GRAPH COORDINATES",
                        required=True)
    parser.add_argument("-g", "--graph", help="AA-formatted graph file", required=True)
    parser.add_argument("--cycle_ids", help="Cycle file elements to use", nargs='+', type=int)
    parser.add_argument("--amp_thresh", help="CN threshold to be considered amplified", type=float, default=10.0)

    args = parser.parse_args()
    elems_to_get = set()
    if args.cycle_ids:
        elems_to_get = set(args.cycle_ids)

    # read the graph
    # retrieve a set of used junctions
    # retrieve a set of segment IDs (mapped to length and CN)
    bidirectional_edge_dict, seg_end_pos_d, segToCN, segToSize = parse_BPG(args.graph)
    nsegs, amp_glen, amp_nsegs = graph_stats(segToSize, segToCN, cutoff=args.amp_thresh)
    seqIDD = id_to_loc(args.graph)

    # take a cycles file and a list of cycle IDs
    # retrieve the cycles as an oriented segment list
    recons_paths = read_cycles(args.cycles, elems_to_get)

    print("Path lengths: ")
    for ind, cp in enumerate(recons_paths):
        path = cp.rsplit(",")
        if path[0] == "0+":
            path = path[1:-1]

        sl = sum([segToSize[int(x[:-1])] for x in path])
        print(str(ind+1) + ": " + str(sl))

    # compute statistics about the amplification of segments
    num_gsegs_exp, prop_gsegs_exp, num_bp_exp, prop_bp_exp = amp_stats(recons_paths, segToSize, amp_glen, amp_nsegs,
                                                                       segToCN, nsegs-1, cutoff=args.amp_thresh)

    # compute statisics about the edges used in AR vs AA
    tot_amp_AA_bps, tot_amp_AR_bps, bp_amp_frac, unused_AA_bps = compare_breakpoints(
        recons_paths, bidirectional_edge_dict, seg_end_pos_d, seqIDD, segToCN, cutoff=args.amp_thresh)

    print("\nComparison of paths " + str(args.cycle_ids) + " between AA and AR\n")
    print("Amplified content used in AA graph vs AR reconstructions")
    print("Total AA amped content: " + str(amp_glen) + " basepairs over " + str(amp_nsegs) + " segments")
    print("AA and AR overlap (bases & proportion) " + str(num_bp_exp) + " (" + str(prop_bp_exp) + ") of AA amped material.")
    print("AA and AR overlap (segments & proportion) " + str(num_gsegs_exp) + " (" + str(prop_gsegs_exp) + ") of AA amped material.\n")

    print("Breakpoints used in AA graph vs AR reconstructions")
    print("Total AA amped breakpoints: " + str(tot_amp_AA_bps))
    print("Total AR amped breakpoints also in AA: " + str(tot_amp_AR_bps) + " (" + str(bp_amp_frac) + ")\n")

    print("Listing of unused amplified breakpoints:")
    for x in unused_AA_bps:
        print(x)