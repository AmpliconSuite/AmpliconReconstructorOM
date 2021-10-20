#!/usr/bin/env python2

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import os
import sys
import argparse
import subprocess

min_map_len = 10  # path must have at least 10 labels for reliable alignment


def run_SegAligner(contig_file, ref_file, arg_list, pre, cnum):
    print("\nCalling SegAligner on " + ref_file)
    cmd_vect = [os.environ['SA_SRC'] + "/SegAligner", ref_file, contig_file] + arg_list
    with open(pre + "_SA_run.out", 'w') as saout:
        print(" ".join(cmd_vect))
        subprocess.call(cmd_vect, stdout=saout)

    aln_files = []
    print("Alignment files generated: ")
    d = "/".join(pre.rsplit("/")[:-1])
    if not d:
        d = "."

    p = pre.rsplit("/")[-1]
    for f in os.listdir(d):
        if f.startswith(p + "_" + str(cnum)) and f.endswith("_aln.txt") and not "rehead" in f:
            ff = d + "/" + f
            print(ff)
            aln_files.append(ff)

    print("")
    return aln_files


def postprocess_alignments(aln_files, cpath, graph_cmap_file):
    print("\nRe-heading and re-indexing alignments based on path")
    for f in aln_files:
        # rewrite header to contain the full path
        of = f.replace("_aln.txt", "_rehead_aln.txt")
        with open(f) as infile, open(of, 'w') as outfile:
            outfile.write(next(infile))
            l = next(infile)
            fields = l.rsplit("\t")
            fields[0] = "#" + cpath
            outfile.write("\t".join(fields))
            for l in infile:
                outfile.write(l)

        # reindex the alignment path based on the segments in the path
        cmd = "{}/scripts/inverse_aln_cycle_ids.py {} {}".format(os.environ["AR_SRC"], graph_cmap_file, of)
        subprocess.call(cmd, shell=True)
        oname = os.path.splitext(os.path.basename(of))[0]
        oname += "_remapped_seg_aln.txt"
        print(oname)
        subprocess.call("rm " + of, shell=True)


def create_path_cmap(cycles_file, gfile, REF, pre, enzyme):
    print("\nConverting cycles to fasta")
    fa_name = ""
    if REF == "GRCh38" or REF == "hg38":
        fa_name = "hg38full.fa"
    elif REF == "hg19":
        fa_name = "hg19full.fa"
    elif REF == "GRCh37":
        fa_name = "human_g1k_v37.fasta"
    else:
        print("ERROR: Could not locate fasta for " + REF)
        sys.exit(1)

    fa_path = "/".join([os.environ['AA_DATA_REPO'], REF, fa_name])
    sp = pre.rstrip("_cycles")
    cmd = os.environ['AR_SRC'] + "/cycles_file_to_fasta.py -c {} -r {} -o {}".format(cycles_file, fa_path, sp)
    print(cmd)
    subprocess.call(cmd, shell=True)

    print("\nCreating graph cmap")
    graph_pre = sp + "_graph"
    cmd = "python2 {}/generate_cmap.py -r {} -e {} -o {} -g {}".format(os.environ['AR_SRC'], fa_path, enzyme, graph_pre,
                                                                       gfile)
    print(cmd)
    subprocess.call(cmd, shell=True)

    print("\nCreating cycles cmap")
    cycles_fa = sp + "_cycles.fasta"
    cmd = "python2 {}/generate_cmap.py -r {} -e {} -o {} --makeRef".format(os.environ['AR_SRC'], cycles_fa, enzyme, pre)
    print(cmd)
    subprocess.call(cmd, shell=True)

    return pre + "_" + enzyme.upper() + ".cmap", graph_pre + "_" + enzyme.upper() + ".cmap"


def get_cycles_from_cnums(cycles_file):
    cycd = {}
    with open(cycles_file) as infile:
        for line in infile:
            if "Cycle=" in line:
                fields = line.rstrip().rsplit(";")
                lineD = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                cnum = int(lineD["Cycle"])
                segs = lineD["Segments"]
                cycd[cnum] = segs

    return cycd


def extract_cmap_by_id(combined_cmap, cnum):
    kf = combined_cmap.replace(".cmap", "_key.txt")
    cmap_id = None
    with open(kf) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.strip().rsplit("\t")
            if fields[1].endswith("_" + str(cnum)):
                cmap_id = fields[0]
                break

    if cmap_id is None:
        print("ERROR: CMAP ID for cycle " + str(cnum) + " not found in key file")
        sys.exit(1)

    cmd = "{}/scripts/extract_cmap.py -i {} -l {}".format(os.environ["AR_SRC"], combined_cmap, cmap_id)
    print(cmd)
    subprocess.call(cmd, shell=True)

    return combined_cmap.rstrip(".cmap") + "_" + cmap_id + ".cmap"


# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Align an AA-formatted path/cycle to a set of OM contigs.")
    parser.add_argument("--contigs", help="contigs cmap file", required=True)
    parser.add_argument("--cycles_file", help="AA-formatted cycles file", required=True)
    parser.add_argument("--graph_file", help="AA-formatted graph file", required=True)
    parser.add_argument("-r", "--ref", help="Reference genome build. Must be same as AA graph",
                        choices=["hg19", "GRCh38", "GRCh37"], required=True)
    parser.add_argument("--cycle_nums", help="IDs of the entry to use from the cycle file (default all)", nargs='+',
                        type=int)
    parser.add_argument("-e", "--enzyme", help="labeling enzyme", choices=["BspQI", "DLE1"], default="DLE1")
    parser.add_argument("-o", "--output_prefix", type=str, help="output alignment filename prefix. Default is cycles "
                                                                "file prefix")
    parser.add_argument("-t", "--threads", help="number of threads to use (default 1)", type=int, default=1)
    parser.add_argument("-i", "--instrument", help="Saphyr or Irys data (default Saphyr)", choices=["Irys", "Saphyr"],
                        default="Saphyr")
    parser.add_argument("--no_clear_alignments", help="Do not clear alignment files from prior runs",
                        action='store_true', default=False)

    args = parser.parse_args()

    gen = "1" if args.instrument == "Irys" else "2"

    # output pathing
    if not args.output_prefix:
        args.output_prefix = os.path.splitext(os.path.basename(args.cycles_file))[0]

    dir_path = "/".join(args.output_prefix.rsplit("/")[:-1])
    try:
        if dir_path and not os.path.exists(dir_path):
            os.makedirs(dir_path)

    except OSError as error:
        print("Directory '%s' can not be created" % dir_path)
        sys.exit(1)

    if not args.no_clear_alignments:
        print("\nRemoving old alignment files from directory.")
        d = dir_path
        if not dir_path:
            d = "."

        p = args.output_prefix.rsplit("/")[-1]
        for f in os.listdir(d):
            if f.startswith(p) and f.endswith("_aln.txt"):
                subprocess.call("rm " + d + "/" + f, shell=True)

    cycd = get_cycles_from_cnums(args.cycles_file)
    if not args.cycle_nums:
        args.cycle_nums = list(cycd.keys())

    # convert cycles to cmap
    combined_cmap_file, graph_cmap_file = create_path_cmap(args.cycles_file, args.graph_file, args.ref,
                                                           args.output_prefix, args.enzyme)

    for cid in sorted(args.cycle_nums):
        individual_cmap = extract_cmap_by_id(combined_cmap_file, cid)
        arg_list = ["-nthreads=" + str(args.threads), "-min_labs=" + str(min_map_len), "-prefix=" + args.output_prefix +
                    "_" + str(cid), "-gen=" + gen]
        aln_files = run_SegAligner(args.contigs, individual_cmap, arg_list, args.output_prefix, cid)
        postprocess_alignments(aln_files, cycd[cid], graph_cmap_file)

    print("\nAll steps completed.")
