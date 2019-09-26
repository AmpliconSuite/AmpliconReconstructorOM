#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import os
import sys
import time
import json
import yaml
import logging
import argparse
import subprocess

try:
	SA_SRC = os.environ['SA_SRC']
	AR_SRC = os.environ['AR_SRC']

except KeyError:
	sys.stderr.write("SA_SRC or AR_SRC bash variable not found. AmpliconReconstructor may not be properly installed.\n")
	sys.exit(1)

def run_ARAD(segs,contigs,graph,enzyme,inst,min_map_len,sname,outdir):
	cmd = "python {}/ARAlignDetect.py -s {} -c {} -g {} -t {} -e {} --min_map_len {} -i {} -o {} -d {}".format(AR_SRC,segs,contigs,graph,nthreads,enzyme,min_map_len,inst,sname,outdir)
	logging.info("CMD:")
	logging.info(cmd)
	subprocess.call(cmd + " > " + outdir + "AS_stdout.txt", shell=True)

def run_OMPF(outdir,segs,contigs,graph,aln_dir,sname,inst,noImpute=False):
	ncs = ""
	if noImpute: ncs+="--noImpute "
	cmd = "python {}/OMPathFinder.py {}--adir {} -c {} -s {} -g {} --outdir {} --prefix {} -i {}".format(AR_SRC,ncs,aln_dir,contigs,segs,graph,outdir,sname,inst)
	subprocess.call(cmd, shell=True)

# def run_visualization(cycles_file,cycle,contigs,segs,graph,aln_file,sname):
# 	# cmd = "python ~/bionano/CycleViz/CycleViz_dev.py --om_alignments --cycles_file {} --cycle {} -c {} -s {} -g {} -i {} --sname {} --label_segs --gene_subset_file {}".format(cycles_file,cycle,contigs,segs,graph,aln_file,sname,os.environ['AR_SRC']+"/Bushman_group_allOnco_May2018.tsv")
# 	subprocess.call(cmd,shell=True)

###"MAIN"###
parser = argparse.ArgumentParser(description="AmpliconReconstructorOM. Wraps methods for alignment and scaffolding.")
parser.add_argument("-i","--yaml_file",type=str, help="Path YAML samples file in AR-OM format")
parser.add_argument("-s",nargs='+', help="Space separated list of samples from YAML dict. If not supplied, all samples are run")
parser.add_argument("--outdir",type=str, help="Directory to store results and associated sub-directories")
parser.add_argument("--run_name",type=str, help="Optional name for this set of runs (will create a directory in outdir of this name)")
parser.add_argument("--nthreads",default=24, type=int, help="number of threads to use")
parser.add_argument("--noAlign",default=False, action='store_true', help="skip alignment step (assume re-using old alignments)")
parser.add_argument("--noImpute",default=False, action='store_true', help="Do not perform path imputation")
parser.add_argument("--noConnect",default=False, action='store_true', help="Do not perform scaffold-linking")
parser.add_argument("--noViz",default=False, action='store_true', help="skip visualizations step (assume re-using old alignments)")
parser.add_argument("--plot_scores",default=False,help="Save plots of the distributions of segment scores")
parser.add_argument("--no_tip_aln",default=False,action='store_true',help="Disable tip alignment step")
parser.add_argument("--no_ref_search",default=False,action='store_true',help="Do not search unaligned regions against reference genome")
parser.add_argument("--no_clear",default=False,action='store_true',help="Do not remove files from previous runs..")
args = parser.parse_args()

if not args.output_directory:
	args.output_directory = os.getcwd()

if not args.outdir.endswith("/"): args.outdir+="/" 
if not args.run_name.endswith("/"): args.run_name+="/" 
run_path = args.outdir + args.run_name
os.mkdir(run_path) if not os.path.exists(run_path)
logging.basicConfig(filename=run_path + "run.log",level=logging.INFO)

noImpute = args.noImpute
noConnect = args.noConnect
nthreads = str(args.nthreads)

logging.info("CMD:")
logging.info(" ".join(sys.argv))
logging.info("noImpute is " + str(noImpute))
logging.info("noConnect is " + str(noImpute))
logging.info("nthreads: " + str(nthreads))

with open(args.yaml_file) as f:
	sample_data = yaml.safe_load(f)

#read samples and figure out what to run (all if no subset specified)
samples_to_run = args.samples if args.samples else [str(x) for x in sample_dict.keys()]

if not args.samples:
	print("running on all samples in " + args.yaml_file)

logging.debug(str(samples_to_run))

#Run the specificed samples
for i in samples_to_run:
	if i not in sample_data:
		logging.warning(i + " not found in YAML file, skipping")
		continue

	logging.info("Processing " + i)

	#ensure the YAML data is correct
	sample_dict = sample_data[i]
	try:
		sample_path = sample_dict["path"]
		seg_path = sample_path + sample_dict["cmap"]
		contig_path = sample_path + sample_dict["contigs"]
		graph_path = sample_path + sample_dict["graph"]
		inst = sample_dict["instrument"]
		enzyme = sample_dict["enzyme"]
		min_map_len = sample_dict["min_map_len"]

	except KeyError:
		em = "YAML file does not contain all required properties for " + i + ", skipping."
		sys.stderr.write(em + "\n")
		logging.warning(em)
		continue

	if any([sample_path,seg_path,contig_path,graph_path,inst,enzyme,min_map_len]) == None:
		em = "Empty property in YAML file for " + i + ", skipping."
		sys.stderr.write(em + "\n")
		logging.warning(em)
		continue

	#make output directories
	alignments_dir = run_path + i + "/alignments/" 
	reconstruction_dir = run_path + i + "/reconstructions/" 
	visualizations_dir = run_path + i + "/visualizations/" 
	logging.info("Making directories")
	for curr_dir in [alignments_dir,reconstruction_dir,visualizations_dir]:
		if not path.exists(curr_dir): os.mkdir(curr_dir) 
		if not args.no_clear:
			logging.info("Clearing old results")
			subprocess.call("rm " + curr_dir + "*",shell=True)


	#remove old "includes detected file"
	idg_basename = os.path.splitext(run_path + i + "/" + sample_dict["graph"])[0]
	#idgf stands for includes_detected graph file
	idgf = graph_basename + "_includes_detected.txt" 
	if os.path.exists(idgf):
		if not args.no_clear:
			logging.warning("Removing old *includes_detected.txt graph file: " + idgf)
			call("rm " + idgf,shell=True)

	start_time = time.time()
	if not args.noAlign:
		run_ARAD(samp_aln_dir, segs_path, contig_path, graph_path, enzyme, inst, i, run_path + i + "/")
	else:
		logging.info("Skipped alignment stage.")
	
	e_time1 = int(time.time() - start_time)
	logging.info("finished alignment stage for " + i + " in " + str(e_time1) + " seconds\nPathfinding")

	#check if output has includes_detected graph file.
	if os.path.exists(idgf):
		idsf = run_path + i + "/" + os.path.splitext(os.path.basename(seg_path))[0] + "_includes_detected.cmap"

	# print segs_path,graph_path
	graph_path,segs_path = idgf,idsf
	logging.info("Using _includes_detected files:\n" + idgf + "\n" + idsf)

	#UPDATE
	run_OMPF(samp_recons_dir, segs_path, contig_path, graph_path, samp_aln_dir + "alignments", i, inst, noImpute)
	e_time2 = (time.time() - start_time) - e_time1
	logging.info("finished pathfinding stage for " + i + " in " + str(e_time2) + " seconds\n")

	# cycles_file = reconstruction_dir + i + "/" + i + "_paths_cycles.txt"
	# cycle = "1" #Need to be able to modify
	# aln_file = reconstruction_dir + i + "/" + i + "_path_" + cycle + "_aln.txt"
	# viz_sname = visualizations_dir + i
	# run_visualization(cycles_file,cycle,contig_path,segs_path,graph_path,aln_file,viz_sname)

	# print("Finished sample " + i)

sys.exit()
