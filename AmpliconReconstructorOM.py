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
import datetime
import subprocess

try:
	SA_SRC = os.environ['SA_SRC']
	AR_SRC = os.environ['AR_SRC']

except KeyError:
	sys.stderr.write("SA_SRC or AR_SRC bash variable not found. AmpliconReconstructor may not be properly installed.\n")
	sys.exit(1)

def run_ARAD(segs,contigs,graph,enzyme,min_map_len,inst,sname,outdir,optionalFlagString=""):
	# psc=""
	# if plot_scores:
	# 	psc+="--plot_scores "
	cmd = "python {}/ARAlignDetect.py {}-s {} -c {} -g {} -t {} -e {} --min_map_len {} -i {} -o {} -d {}".format(AR_SRC,optionalFlagString,segs,contigs,graph,nthreads,enzyme,min_map_len,inst,sname,outdir)
	logging.info("ARAD CMD:")
	logging.info(cmd)
	subprocess.call(cmd, shell=True)

def run_OMPF(outdir,segs,contigs,graph,aln_dir,sname,inst,optionalFlagString = ""):
	cmd = "python {}/OMPathFinder.py {}--adir {} -c {} -s {} -g {} --outdir {} --prefix {} -i {}".format(AR_SRC,optionalFlagString,aln_dir,contigs,segs,graph,outdir,sname,inst)
	logging.info("OMPF CMD:")
	logging.info(cmd)
	subprocess.call(cmd,shell=True)

# def run_visualization(cycles_file,cycle,contigs,segs,graph,aln_file,sname):
# 	# cmd = "python ~/bionano/CycleViz/CycleViz_dev.py --om_alignments --cycles_file {} --cycle {} -c {} -s {} -g {} -i {} --sname {} --label_segs --gene_subset_file {}".format(cycles_file,cycle,contigs,segs,graph,aln_file,sname,os.environ['AR_SRC']+"/Bushman_group_allOnco_May2018.tsv")
# 	subprocess.call(cmd,shell=True)

###"MAIN"###
parser = argparse.ArgumentParser(description="AmpliconReconstructorOM. Wraps methods for alignment and scaffolding.")
parser.add_argument("-i","--yaml_file",type=str, help="Path YAML samples file in AR-OM format",required=True)
parser.add_argument("-s","--samples",nargs='+', help="If only running a subset, -s followed by space separated list of samples names from YAML \
					dict. If not supplied, all samples are run")
parser.add_argument("--outdir",type=str, help="Directory to store results and associated sub-directories",required=True)
parser.add_argument("--run_name",type=str, default="", help="Optional name for this set of runs (will create a directory in outdir of this name)")
parser.add_argument("--nthreads",default=24, type=int, help="number of threads to use")
parser.add_argument("--noAlign",default=False, action='store_true', help="skip alignment step (assume re-using old alignments)")
parser.add_argument("--noImpute",default=False, action='store_true', help="Do not perform path imputation")
parser.add_argument("--noConnect",default=False, action='store_true', help="Do not perform scaffold-linking")
parser.add_argument("--noViz",default=False, action='store_true', help="skip visualizations step (assume re-using old alignments)")
parser.add_argument("--plot_scores",default=False,action='store_true',help="Save plots of the distributions of segment scores")
parser.add_argument("--no_tip_aln",default=False,action='store_true',help="Disable tip alignment step")
parser.add_argument("--no_ref_search",default=False,action='store_true',help="Do not search unaligned regions against reference genome")
parser.add_argument("--no_clear",default=False,action='store_true',help="Do not remove files from previous runs..")
args = parser.parse_args()

if not args.outdir:
	args.outdir = os.getcwd()

if not args.outdir.endswith("/"): args.outdir+="/" 
if not args.run_name.endswith("/"): args.run_name+="/" 
run_path = args.outdir + args.run_name
if not os.path.exists(run_path): os.mkdir(run_path) 
logging.basicConfig(filename=run_path + "run.log",level=logging.INFO)
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

#read samples and figure out what to run (all if no subset specified)
samples_to_run = args.samples if args.samples else [str(x) for x in sample_dict.keys()]

if not args.samples: print("running on all samples in " + args.yaml_file)

logging.debug(str(samples_to_run))

#Run the specificed samples
for i in samples_to_run:
	if i not in sample_data:
		logging.error(i + " not found in YAML file, skipping")
		continue

	logging.info("\nProcessing " + i)
	print("\nProcessing " + i)

	#ensure the YAML data is correct
	sample_dict = sample_data[i]
	try:
		sample_path = sample_dict["path"]
		segs_path = sample_path + sample_dict["cmap"]
		contig_path = sample_path + sample_dict["contigs"]
		graph_path = sample_path + sample_dict["graph"]
		inst = sample_dict["instrument"]
		enzyme = sample_dict["enzyme"]
		min_map_len = sample_dict["min_map_len"]

	except KeyError:
		em = "YAML file does not contain all required properties for " + i + ", skipping."
		sys.stderr.write(em + "\n")
		logging.error(em)
		continue

	if any([sample_path,segs_path,contig_path,graph_path,inst,enzyme,min_map_len]) == None:
		em = "Empty property in YAML file for " + i + ", skipping."
		sys.stderr.write(em + "\n")
		logging.error(em)
		continue

	#make output directories
	rpi = run_path + i + "/"
	alignments_dir = rpi + "alignments/" 
	reconstruction_dir = rpi + "reconstructions/" 
	visualizations_dir = rpi + "visualizations/" 
	logging.info("Making directories")
	for curr_dir in [rpi,alignments_dir,reconstruction_dir,visualizations_dir]:
		if not os.path.exists(curr_dir): os.mkdir(curr_dir) 
		if not args.no_clear:
			logging.info("Clearing old results")
			subprocess.call("rm " + curr_dir + "* 2>/dev/null",shell=True)


	#remove old "includes detected file"
	#idgf stands for includes_detected graph file
	idgf = os.path.splitext(rpi + sample_dict["graph"])[0] + "_includes_detected.txt"
	if os.path.exists(idgf):
		if not args.no_clear:
			logging.warning("Removing old *includes_detected.txt graph file: " + idgf)
			call("rm " + idgf,shell=True)

	#check if we're converting an xmap:
	use_xmap = False
	if "xmap" in sample_dict:
		if sample_dict["xmap"]: use_xmap = True

	start_time = time.time()
	optionalFlagString = ""
	if not args.noAlign:
		if args.plot_scores:
			optionalFlagString+="--plot_scores "
		run_ARAD(segs_path, contig_path, graph_path, enzyme, min_map_len, inst, i, rpi, optionalFlagString)

	elif use_xmap: #elif check xmap status
		optionalFlagString+="--xmap " + sample_dict["xmap"] + " "
		run_ARAD(segs_path, contig_path, graph_path, enzyme, min_map_len, inst, i, rpi, optionalFlagString)

	else:
		logging.info("Skipped alignment stage.")
	
	e_time1 = int(time.time() - start_time)
	logging.info("finished alignment stage for " + i + " in " + str(e_time1) + " seconds\nPathfinding")

	#check if output has includes_detected graph file.
	if os.path.exists(idgf):
		idsf = rpi + os.path.splitext(os.path.basename(segs_path))[0] + "_includes_detected.cmap"
		graph_path,segs_path = idgf,idsf
		logging.info("Using _includes_detected files:\n" + idgf + "\n" + idsf)

	print "Reconstructing amplicon " + i
	optionalFlagString = ""
	if args.noImpute:
		optionalFlagString+="--noImpute "

	if args.noConnect:
		optionalFlagString+="--noConnect "

	run_OMPF(reconstruction_dir, segs_path, contig_path, graph_path, alignments_dir, i, inst, optionalFlagString)
	e_time2 = (time.time() - start_time) - e_time1
	logging.info("finished pathfinding stage for " + i + " in " + str(e_time2) + " seconds\n")

	cycles_file = reconstruction_dir + i + "_paths_cycles.txt"
	#determine number of paths

	#output the length of each path and each scaffold 

	#for path in paths
		#make linear
		#make cyclic
		#get associated _aln.txt file
		

	# cycle = "1" #Need to be able to modify
	# aln_file = reconstruction_dir + i + "/" + i + "_path_" + cycle + "_aln.txt"
	# viz_sname = visualizations_dir + i
	# run_visualization(cycles_file,cycle,contig_path,segs_path,graph_path,aln_file,viz_sname)

	# print("Finished sample " + i)

logging.info("Finished")
logging.info(str(datetime.datetime.now()) + "\n")
logging.shutdown()
sys.exit()
