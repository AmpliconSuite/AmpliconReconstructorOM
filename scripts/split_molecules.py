import sys
import os
import copy
import logging
import argparse
import numpy as np
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from bionanoUtil import *
from readclust import *

clustDelta = 50000
min_clust_size = 3
aln_cut_percentile = 50 #percentile of cut

def read_graph(graphf):
    gseqs = defaultdict(IntervalTree)
    deList = []

    with open(graphf) as infile:
        for line in infile:
            if line.startswith("sequence"):
                fields = line.rstrip().rsplit()
                chrom = fields[1].split(":")[0]
                p1 = int(fields[1].split(":")[1][:-1])
                p2 = int(fields[2].split(":")[1][:-1])
                gseqs[chrom][p1:p2] = float(fields[3])

            if line.startswith("discordant"):
                fields = line.rstrip().rsplit()
                lbp,rbp = fields[1].split("->")
                lchrom,lpd = lbp.rsplit(":")
                rchrom,rpd = rbp.rsplit(":")
                lpos,ldir = int(lpd[:-1]),lpd[-1]
                rpos,rdir = int(rpd[:-1]),rpd[-1]
                rSupp = int(fields[3])
                # isFoldback = (ldir == rdir)

                r1 = dummy_read(lchrom, lpos, ldir == "-")
                r2 = dummy_read(rchrom, rpos, rdir == "-")
                # sr1,sr2 = sorted([r1,r2],key=lambda x: (x.reference_name,x.reference_end))

                curr_clust = pe_read_clust(r1,r2,clustDelta)
                print(str(curr_clust.clust_to_bedpe()))
                curr_clust.size = rSupp
                deList.append(curr_clust)

    return gseqs,deList


def read_excludedRegions(exc_file,ref):
    excIT = defaultdict(IntervalTree)
    with open(exc_file) as infile:
        for line in infile:
            fields = line.rstrip().rsplit("\t")
            fields[1], fields[2] = int(fields[1]), int(fields[2])
            if ref == "GRCh37" and fields[0].startswith("chr"):
                fields[0] = fields[0][3:]

            excIT[fields[0]].add(Interval(fields[1], fields[2]))

    return excIT


def clustIsExcludeable(excIT, clust):
    ref1, ref2 = clust.r_IDs[0], clust.r_IDs[0]
    p1, p2 = clust.centroid
    return excIT[ref1][p1] or excIT[ref2][p2]


# match seqs on length
def bionanoID_to_chrom(ref_lens, ref_key):
    id_chrom_map = {}
    len_to_key = {v[1]: v[0] for v in ref_key.values()}
    for k, s in ref_lens.iteritems():
        try:
            currChrom = len_to_key[s]
            id_chrom_map[k] = currChrom

        except KeyError:
            sys.stderr.write("Warning: entry " + k + " with length " + str(s) + " not found in ref key\n")
            continue

    return id_chrom_map


# check if segment overlaps bed
def segment_bed_match(chrom, region, bed_dict):
    if chrom not in bed_dict:
        return None

    return bed_dict[chrom][region[0]:region[1]]


def xmap_contains_multimapping(xmap_d):
    seen_qrys = set()
    for k, v in xmap_d.iteritems():
        qry_id = v["QryContigID"]
        if qry_id in seen_qrys:
            return True

        else:
            seen_qrys.add(qry_id)

    return False


# methods for bp extn
def get_mols_with_bed_aln(xmap_d, bed_dict, id_chrom_map):
    mol_to_aln_list = defaultdict(list)
    rel_mols = set()
    for _, aln_d in xmap_d.iteritems():
        m_id = aln_d["QryContigID"]
        mol_to_aln_list[m_id].append(aln_d)
        curr_chrom = id_chrom_map[aln_d["RefContigID"]]
        ar_s = aln_d["RefStartPos"]
        ar_e = aln_d["RefEndPos"]

        hit_region = segment_bed_match(curr_chrom, (ar_s, ar_e), bed_dict)
        if hit_region:
            rel_mols.add(m_id)

    rel_mol_sorted_alns = {x: sorted(mol_to_aln_list[x], key=lambda y: min(y["QryStartPos"], y["QryEndPos"])) for x in
                           rel_mols}


    return rel_mol_sorted_alns


def get_aln_norm_score(aln_d):
    qStart = aln_d["Alignment"][0][1]
    qEnd = aln_d["Alignment"][-1][1]
    return aln_d["Confidence"]/abs(qEnd-qStart)


def aln_percentile_score(xmap_d,percentile):
    norm_scores = []
    for _, aln_d in xmap_d.iteritems():
        qStart = aln_d["Alignment"][0][1]
        qEnd = aln_d["Alignment"][-1][1]
        norm_scores.append(get_aln_norm_score(aln_d))

    sn_scores = sorted(norm_scores)
    #return value of list at percentile
    aps = np.percentile(sn_scores,percentile)
    logging.debug("Percentile value" + str(aps))
    return aps


#return inSegs,inGraph
def pe_read_in_graph(r1, r2, gseqs, deList):
    chrom1, s1, e1 = r1.reference_name, r1.reference_start, r1.reference_end
    chrom2, s2, e2 = r2.reference_name, r2.reference_start, r2.reference_end
    relSegInts1 = gseqs[chrom1][s1:e1]
    relSegInts2 = gseqs[chrom2][s2:e2]
    inSegs = int(len(relSegInts1) > 0) + int(len(relSegInts2) > 0)
    for gc in deList:
        if gc.rp_has_overlap(r1, r2):
            return inSegs, True

    return inSegs, False


def clust_in_graph(cc, gseqs, deList):
    chrom1, chrom2 = cc.r_IDs
    s1, e1 = cc.centroid
    relSegInts1 = gseqs[chrom1][s1]
    relSegInts2 = gseqs[chrom2][e1]
    inSegs = int(len(relSegInts1) > 0) + int(len(relSegInts2) > 0)
    for gc in deList:
        if gc.clust_has_overlap(cc):
            return inSegs, True

    return inSegs, False


def filter_rel_mol_alns(rel_mol_sorted_alns,molLenD,aln_cut,excIT,id_chrom_map):
    #check if it has an overarching high-confidence alignment
    #define high confidence as 50% or greater average score
    filt_mol_sorted_alns = {}
    for m_id, aln_list in rel_mol_sorted_alns.iteritems():
        if len(aln_list) > 1: #ONLY EXAMINE MOLS WITH > 1 ALN
            has_overarching_aln = False
            for xd in aln_list:
                mol_aln_len = abs(xd["QryEndPos"] - xd["QryStartPos"])
                mol_percent_aln = mol_aln_len/molLenD[m_id]
                aln_norm_score = get_aln_norm_score(xd)
                curr_chrom = id_chrom_map[xd["RefContigID"]]
                if excIT[curr_chrom][(xd["RefStartPos"],xd["RefEndPos"])]:
                    continue

                if mol_percent_aln > 0.8 and aln_norm_score > aln_cut:
                    has_overarching_aln = True
                    break

            if not has_overarching_aln:
                filt_mol_sorted_alns[m_id] = aln_list

    return filt_mol_sorted_alns


#turn molecules into "reads"
def mol_alns_to_read(filt_mol_sorted_alns,id_chrom_map):
    #turn the molecule alignments into a list of dummy reads
    dummy_reads = []
    for m_id, aln_list in filt_mol_sorted_alns.iteritems():
        dummy_reads.append([])
        for xd in aln_list:
            curr_chrom = id_chrom_map[xd["RefContigID"]]
            curr_dummy_read = dummy_read(curr_chrom, xd["RefStartPos"], xd["Orientation"], m_id)
            curr_dummy_read.qstart,curr_dummy_read.qend = xd["QryStartPos"],xd["QryEndPos"]
            curr_dummy_read.reference_end = xd["RefEndPos"]
            curr_dummy_read.is_read2 = False
            curr_dummy_read.is_supplemental = True
            curr_dummy_read.mapping_quality = get_aln_norm_score(xd)
            dummy_reads[-1].append(curr_dummy_read)

    return dummy_reads


#create a "sortedv" which is generated by pairwise matchings of the dummy reads
def mol_reads_to_sorted_read_pairs(unpaired_dummy_reads):
    unsorted_pairs = defaultdict(list)
    for rl in unpaired_dummy_reads:
        if len(rl) > 1: #only keep the entries with more than one read
            for i in range(len(rl)):
                for j in range(i,len(rl)):
                    i_true_end = max(rl[i].qstart,rl[i].qend)
                    j_true_start = min(rl[j].qstart,rl[j].qend)

                    if j_true_start + 10000 < i_true_end: #too much overlap
                        continue

                    vp = (rl[i],rl[j])
                    sortedv = tuple(sorted(vp, key=lambda x: (x.reference_name, x.reference_end)))
                    refpair = (sortedv[0].reference_name,sortedv[1].reference_name)
                    unsorted_pairs[refpair].append(sortedv)

    final_pairs = dict()
    for k,l in unsorted_pairs.iteritems():
        sorted_reads = sorted(l,key=lambda x: (x[0].reference_end, x[1].reference_start))
        final_pairs[k] = sorted_reads

    return final_pairs


def cluster_discordant_reads(sr_dict,excIT):
    clusts = defaultdict(list)
    print(len(sr_dict))
    for cp,sr in sr_dict.items():
        print(cp,len(sr))
        curr_clusts = []
        curr_clust = pe_read_clust(sr[0][0],sr[0][1],clustDelta)
        curr_clusts.append(curr_clust)
        for r1,r2 in sr[1:]:
            found = False

            for cc in curr_clusts:
                if r1.reference_end - cc.centroid[0] < 2*clustDelta:
                    if cc.rp_has_overlap(r1,r2):
                        cc.add_pair_to_clust(r1,r2)
                        found = True
                        break

            if not found:
                curr_clusts.append(pe_read_clust(r1,r2,clustDelta))

        #fencepost at end
        for cc in curr_clusts:
            if cc.size >= min_clust_size:
                if not clustIsExcludeable(excIT,cc):
                    clusts[cp].append(copy.copy(cc))

    return clusts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Identify bionano breakpoints involving AA graph segments")
    parser.add_argument("--ref_version", help="Reference genome version.", choices=["hg19", "GRCh37", "GRCh38"],
                        required=True)
    parser.add_argument("--ref_fname", help="Alignment reference file", required=True)
    parser.add_argument("--query_fname", help="Alignment query file", required=True)
    parser.add_argument("--xmap_fname",help="Alignment xmap file", required=True)
    parser.add_argument("--AR_ref_key",help="Path to the AR ref '_key.txt' file for same genome version used with "
                                            "'--ref_fname'")
    parser.add_argument("--exclude",help="File of breakpoint excludable regions",required=True)
    parser.add_argument("--AA_graph",help="path to AA graph",required=True)
    parser.add_argument("-o",help="Output filename prefix")
    args = parser.parse_args()

    base = os.path.basename(args.AA_graph)
    basename = os.path.splitext(base)[0]
    if not args.o:
        args.o = basename

    logging.basicConfig(filename=args.o + '_splitmol.log',
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=logging.DEBUG,
                        datefmt='%Y-%m-%d %H:%M:%S')


    graph_seg_dict, deList = read_graph(args.AA_graph)


    logging.info("Getting mol_lens")
    molLenD = get_cmap_lens(args.query_fname)

    logging.info("parsing ref, key & making vects")
    ref_cmaps = parse_cmap(args.ref_fname)

    logging.info("getting excludeable regions")
    excIT = read_excludedRegions(args.exclude, args.ref_version)

    # this should make the key from the other key
    ref_key = parse_keyfile(args.AR_ref_key)
    ref_lens = get_cmap_lens(args.ref_fname)
    id_chrom_map = bionanoID_to_chrom(ref_lens, ref_key)
    ref_vects = vectorize_cmaps(ref_cmaps)

    logging.info("parsing xmap")
    xmap_d = parse_xmap(args.xmap_fname)
    logging.info(str(len(xmap_d)) + " entries")
    logging.info("XMAP is multimapping: " + str(xmap_contains_multimapping(xmap_d)))

    #get threshold for confident align
    aln_cut = aln_percentile_score(xmap_d,aln_cut_percentile)

    #get relevant sorted molecules
    rel_mol_sorted_alns = get_mols_with_bed_aln(xmap_d, graph_seg_dict, id_chrom_map)
    logging.debug("Number of rel mol ids (pre-filter): " + str(len(rel_mol_sorted_alns)))

    #filter the molecules
    filt_mol_sorted_alns = filter_rel_mol_alns(rel_mol_sorted_alns, molLenD, aln_cut, excIT, id_chrom_map)
    logging.debug("Number of rel mol ids (post-filter): " + str(len(filt_mol_sorted_alns)))

    #get an unsorted list of dummy reads
    unpaired_dummy_reads = mol_alns_to_read(filt_mol_sorted_alns, id_chrom_map)

    #pair up the dummy reads into dummy discordant read pairs
    final_pairs = mol_reads_to_sorted_read_pairs(unpaired_dummy_reads)

    #write some raw reads
    logging.info("Writing raw read output")
    with open(args.o + "_mol_raw_discordant.bedpe",'w') as mol_of:
        mol_of.write("Sample\tLeftChr\tLeftEnd\tRightChr\tRightStart\tMolID\tinSegs\tinGraph\n")
        for cp, rp_l in final_pairs.items():
            for r1, r2 in rp_l:
                if r1.query_name != r2.query_name:
                    logging.warning("Read pair had non-matching query_names " + r1.query_name + " " + r2.query_name)
                    continue
                inSegs, inGraph = pe_read_in_graph(r1, r2, graph_seg_dict, deList)
                mol_of.write("\t".join([str(x) for x in [basename, r1.reference_name, r1.reference_end, r2.reference_name,
                                                         r2.reference_start, r1.query_name, inSegs, inGraph]]) + "\n")

    #cluster the dummy discordant read pairs
    sDR_clusts = cluster_discordant_reads(final_pairs, excIT)
    logging.info("Writing read clusts")
    with open(args.o + "_mol_discordant_clusts.txt", 'w') as of1, open(args.o + "_mol_clust_read_info.txt",'w') as of2:
        #iterate over clusts
        of1.write("Sample\tstart_chr\tstart_pos\tend_chr\tend_pos\tNumReads\tinSegs\tinGraph\n")
        for k,l in sDR_clusts.items():
            for cc in l:
                #check if a cluster matches something in the graph file.
                inSegs, inGraph = clust_in_graph(cc, graph_seg_dict, deList)
                #write each clust to file
                of1.write("\t".join([str(x) for x in [basename,] + cc.clust_to_bedpe() + [inSegs,inGraph]]) + "\n")
                of2.write(cc.clust_to_string() + "\n")

    logging.info("Finished\n")
    print("Finished")
    sys.exit()






