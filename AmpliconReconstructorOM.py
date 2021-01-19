#!/usr/bin/env python2

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

__author__ = "Jens Luebeck"
__version__ = "1.01"

import os
import sys
import time
import json
import yaml
import logging
import argparse
import datetime
import subprocess

try:
    SA_SRC = os.environ['SA_SRC']
    AR_SRC = os.environ['AR_SRC']

except KeyError:
    sys.stderr.write("SA_SRC or AR_SRC bash variable not found. AmpliconReconstructor may not be properly installed.\n")
    sys.exit(1)


def run_ARAD(segs, contigs, graph, enzyme, min_map_len, inst, sname, outdir, ref, optionalFlagString=""):
    # psc=""
    # if plot_scores:
    # 	psc+="--plot_scores "
    if min_map_len is None:
        cmd = "python {}/ARAlignDetect.py {}-s {} -c {} -g {} -t {} -e {} -i {} -o {} -d {} -r {}".format(
            AR_SRC, optionalFlagString, segs, contigs, graph, nthreads, enzyme, inst, sname, outdir, ref)

    else:
        cmd = "python {}/ARAlignDetect.py {}-s {} -c {} -g {} -t {} -e {} --min_map_len {} -i {} -o {} -d {} -r {}".format(
            AR_SRC, optionalFlagString, segs, contigs, graph, nthreads, enzyme, min_map_len, inst, sname, outdir, ref)

    logging.info("ARAD CMD:")
    logging.info(cmd)
    subprocess.call(cmd, shell=True)


def run_OMPF(outdir, segs, contigs, graph, aln_dir, sname, inst, optionalFlagString=""):
    cmd = "python {}/OMPathFinder.py {}--adir {} -c {} -s {} -g {} --outdir {} --prefix {} -i {}".format(
        AR_SRC, optionalFlagString, aln_dir, contigs, segs, graph, outdir, sname, inst)

    logging.info("OMPF CMD:")
    logging.info(cmd)
    subprocess.call(cmd, shell=True)


def parse_cycles_file(cycles_file):
    # parse the file(s)
    segs = {}
    paths = {}
    with open(cycles_file) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                fnum = fields[1]
                chrom, p1, p2 = fields[2], int(fields[3]), int(fields[4])
                segs[fnum] = (chrom, p1, p2)

            elif line.startswith("Cycle"):
                cycleFields = line.rstrip().rsplit(";")
                cycleD = {}
                for f in cycleFields:
                    n, v = f.rsplit("=")
                    cycleD[n] = v

                paths[cycleD["Cycle"]] = cycleD["Segments"].rsplit(",")

    return segs, paths


def compute_path_lengths(segs, paths):
    lens = {}
    for k, v in paths.items():
        totlen = 0.0
        for x in v:
            if x[:-1] != "0":
                stup = segs[x[:-1]]
                slen = stup[2] - stup[1]
                totlen += slen

        lens[k] = totlen

    return lens


def run_visualization(CV_path, cycles_file, cycleNum, contigs, segs, graph, aln, sname, label_segs=True,
                      subset_genes=True):
    optionalFlagStringCV = ""
    optionalFlagStringLV = ""
    if label_segs:
        optionalFlagStringCV += "--label_segs numbers "
        optionalFlagStringLV += "--label_segs id "

    if subset_genes:
        gspath = "--gene_subset_file Bushman "
        optionalFlagStringCV += gspath
        optionalFlagStringLV += gspath

    cmd = "python {}/CycleViz.py --om_alignments {}--cycles_file {} --cycle {} -c {} --om_segs {} -g {} -i {} --outname {}".format(
        CV_path, optionalFlagStringCV, cycles_file, cycleNum, contigs, segs, graph, aln, sname)
    subprocess.call(cmd, shell=True)

    cmd = "python {}/LinearViz.py --om_alignments {}--cycles_file {} --path {} -c {} --om_segs {} -g {} -i {} --outname {}".format(
        CV_path, optionalFlagStringLV, cycles_file, cycleNum, contigs, segs, graph, aln, sname)
    subprocess.call(cmd, shell=True)


###"MAIN"###
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="AmpliconReconstructorOM. Wraps methods for alignment and scaffolding.")
    parser.add_argument("-i", "--yaml_file", type=str, help="Path YAML samples file in AR-OM format", required=True)
    parser.add_argument("-s", "--samples", nargs='+', help="If only running a subset, -s followed by space separated "
                                                           "list of samples names from YAML dict. If not supplied, "
                                                           "all samples are run")
    parser.add_argument("--outdir", type=str, help="Directory to store results and create associated sub-directories",
                        required=True)
    parser.add_argument("--run_name", type=str, default="",
                        help="Optional name for this set of runs (will create a directory in outdir of this name)")
    parser.add_argument("--nthreads", default=24, type=int, help="number of threads to use")
    parser.add_argument("--noAlign", default=False, action='store_true',
                        help="skip alignment step (assume re-using old alignments)")
    parser.add_argument("--noImpute", default=False, action='store_true', help="Do not perform path imputation")
    parser.add_argument("--noConnect", default=False, action='store_true', help="Do not perform scaffold-linking")
    parser.add_argument("--noViz", default=False, action='store_true', help="skip visualizations step")
    parser.add_argument("--plot_scores", default=False, action='store_true',
                        help="Save plots of the distributions of segment scores")
    parser.add_argument("--no_tip_aln", default=False, action='store_true', help="Disable tip alignment step")
    parser.add_argument("--no_ref_search", default=False, action='store_true',
                        help="Do not search unaligned regions against reference genome")
    parser.add_argument("--no_clear", default=False, action='store_true',
                        help="Do not remove files from previous runs..")
    parser.add_argument("--visualize_scaffolds", default=False, action='store_true',
                        help="Produce visualizations for the single unconnected scaffolds")
    parser.add_argument("--CV_path", help="Path to the CycleViz source directory")
    args = parser.parse_args()

    if not args.outdir:
        args.outdir = os.getcwd()

    elif not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    if not args.outdir.endswith("/"): args.outdir += "/"
    if not args.run_name.endswith("/"): args.run_name += "/"
    run_path = args.outdir + args.run_name
    if not os.path.exists(run_path): os.mkdir(run_path)
    logging.basicConfig(filename=run_path + "run.log", level=logging.INFO, filemode='w')
    logging.info(str(datetime.datetime.now()))
    logging.info("Starting logging")

    noImpute = args.noImpute
    noConnect = args.noConnect
    nthreads = str(args.nthreads)

    logging.info("CMD:")
    logging.info(" ".join(sys.argv))
    logging.info("noImpute is " + str(noImpute))
    logging.info("noConnect is " + str(noConnect))
    logging.info("nthreads: " + str(nthreads))

    with open(args.yaml_file) as f:
        sample_data = yaml.safe_load(f)

    # read samples and figure out what to run (all if no subset specified)
    samples_to_run = args.samples if args.samples else [str(x) for x in sample_data.keys()]

    if not args.samples: print("running on all samples in " + args.yaml_file)

    logging.debug(str(samples_to_run))

    # Run the specificed samples
    for i in samples_to_run:
        if i not in sample_data:
            logging.error(i + " not found in YAML file, skipping")
            continue

        logging.info("\nProcessing " + i)
        print("\nProcessing " + i)

        # ensure the YAML data is correct
        sample_dict = sample_data[i]
        try:
            sample_path = sample_dict["path"]
            if not sample_path.endswith("/"): sample_path += "/"
            segs_path = sample_path + sample_dict["cmap"]
            contigs_path = sample_path + sample_dict["contigs"]
            graph_path = sample_path + sample_dict["graph"]
            inst = sample_dict["instrument"]
            enzyme = sample_dict["enzyme"]
            min_map_len = sample_dict["min_map_len"]
            ref = sample_dict["reference_build"]
            if ref.lower() == "hg38": ref = "GRCh38"

        except KeyError:
            em = "YAML file does not contain properties for " + i + ". Skipping sample."
            sys.stderr.write(em + "\n")
            logging.error(em)
            continue

        if any([x is None for x in [sample_path, segs_path, contigs_path, graph_path, inst, enzyme, ref]]):
            em = "None-type in required property in YAML file for " + i + ". Skipping sample."
            sys.stderr.write(em + "\n")
            logging.error(em)
            continue

        # make output directories
        rpi = run_path + i + "/"
        alignments_dir = rpi + "alignments/"
        reconstruction_dir = rpi + "reconstructions/"
        visualizations_dir = rpi + "visualizations/"
        logging.info("Making directories")
        for curr_dir in [rpi, alignments_dir, reconstruction_dir, visualizations_dir]:
            if not os.path.exists(curr_dir): os.mkdir(curr_dir)
            if not args.no_clear and not (curr_dir == alignments_dir and args.noAlign) and not curr_dir == rpi:
                logging.info("Clearing old results")
                subprocess.call("rm " + curr_dir + "* 2>/dev/null", shell=True)

        # remove old "includes detected file"
        # idgf stands for includes_detected graph file
        idgf = os.path.splitext(rpi + sample_dict["graph"].rsplit("/")[-1])[0] + "_includes_detected.txt"
        # print idgf,os.path.exists(idgf)
        if os.path.exists(idgf):
            if not args.no_clear:
                print("Removing old *includes_detected.txt graph file: " + idgf)
                subprocess.call("rm " + idgf, shell=True)

        # check if we're converting an xmap:
        use_xmap = False
        if "xmap" in sample_dict:
            if sample_dict["xmap"]: use_xmap = True

        start_time = time.time()

        # RUN ARAD
        optionalFlagString = ""
        if not args.noAlign:
            if args.plot_scores:
                optionalFlagString += "--plot_scores "

            if args.no_ref_search:
                optionalFlagString += "--no_ref_search "

            if args.no_tip_aln:
                optionalFlagString += "--no_tip_aln "

            run_ARAD(segs_path, contigs_path, graph_path, enzyme, min_map_len, inst, i, rpi, ref, optionalFlagString)

        elif use_xmap:  # elif check xmap status
            optionalFlagString += "--xmap " + sample_dict["xmap"] + " "
            run_ARAD(segs_path, contigs_path, graph_path, enzyme, min_map_len, inst, i, rpi, ref, optionalFlagString)

        else:
            logging.info("Skipped alignment stage.")

        e_time1 = int(time.time() - start_time)
        logging.info("finished alignment stage for " + i + " in " + str(e_time1) + " seconds\nPathfinding")

        # check if output has includes_detected graph file.
        if os.path.exists(idgf):
            idsf = rpi + os.path.splitext(os.path.basename(segs_path))[0] + "_includes_detected.cmap"
            graph_path, segs_path = idgf, idsf
            logging.info("Using _includes_detected files:\n" + idgf + "\n" + idsf)

        # RUN OMPF
        print("Reconstructing amplicon " + i)
        optionalFlagString = ""
        if args.noImpute:
            optionalFlagString += "--noImpute "

        if args.noConnect:
            optionalFlagString += "--noConnect "

        run_OMPF(reconstruction_dir, segs_path, contigs_path, graph_path, alignments_dir, i, inst, optionalFlagString)
        e_time2 = (time.time() - start_time) - e_time1
        logging.info("finished pathfinding stage for " + i + " in " + str(e_time2) + " seconds\n")
        logging.info("doing visualizations")

        paths_cycles_file = reconstruction_dir + i + "_paths_cycles.txt"
        scaffolds_file = reconstruction_dir + i + "_scaffold_paths.txt"

        # PATH LENGTH AND VISUALIZATION
        # determine number of paths
        path_alns = [x for x in os.listdir(reconstruction_dir) if i + "_path_" in x]
        pnums = [x.rsplit(i + "_path_")[1].rsplit("_")[0] for x in path_alns]
        pnum_to_alnfile = dict(zip(pnums, path_alns))
        scaffold_alns = [x for x in os.listdir(reconstruction_dir) if i + "_scaffold_path_" in x]
        snums = [x.rsplit(i + "_scaffold_path_")[1].rsplit("_")[0] for x in scaffold_alns]
        snum_to_alnfile = dict(zip(snums, scaffold_alns))

        # output the length of each path and each scaffold
        path_segs, paths = parse_cycles_file(paths_cycles_file)
        path_lens = compute_path_lengths(path_segs, paths)

        scaffold_segs, scaffolds = parse_cycles_file(scaffolds_file)
        scaffold_lens = compute_path_lengths(scaffold_segs, scaffolds)

        with open(reconstruction_dir + i + "_scaffolds_cycles_lengths.txt", 'w') as outfile:
            outfile.write("#Path lengths\n")
            for ki in sorted([int(x) for x in paths.keys()]):
                k = str(ki)
                outfile.write("\t".join([k, ",".join(paths[k]), str(path_lens[k])]) + "\n")

            outfile.write("#Scaffold lengths\n")
            for ki in sorted([int(x) for x in scaffolds.keys()]):
                k = str(ki)
                outfile.write("\t".join([k, ",".join(scaffolds[k]), str(scaffold_lens[k])]) + "\n")

        if not args.CV_path:
            try:
                args.CV_path = os.environ['CV_SRC']
            except KeyError:
                pass

        if not args.noViz and args.CV_path:
            for k in sorted(paths.keys())[:40]:
                aln_file = reconstruction_dir + pnum_to_alnfile[k]
                run_visualization(args.CV_path, paths_cycles_file, k, contigs_path, segs_path, graph_path, aln_file,
                                  visualizations_dir + i)

            logging.info("finished visualization stage for " + i + "\n")

    logging.info("Finished")
    logging.info(str(datetime.datetime.now()) + "\n")
    logging.shutdown()
    sys.exit()
