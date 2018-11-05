import sys
import os
import time
import bisect
import argparse
import threading
import subprocess
import numpy as np
from scipy.stats import linregress
from collections import defaultdict
from bionanoUtil import *

#development version

#coordinates the pipeline

unaligned_label_cutoff = 8
unaligned_size_cutoff = 75000

#break the cmap into smaller files
def chunk_cmap(cmap_file,nthreads,min_lab=0):
    cmap_header = []
    cmap_data = defaultdict(list)
    cmap_chunked_fnames = []
    with open(cmap_file) as infile:
        for line in infile:
            if line.startswith("#"):
                cmap_header.append(line)

            else:
                cmap_id = int(line.rsplit("\t")[0])
                cmap_data[cmap_id].append(line)

    keys_to_del = set()
    for i in cmap_data:
        if len(cmap_data[i]) < min_lab:
            keys_to_del.add(i)

    for i in keys_to_del:
        cmap_data.pop(i,None)

    chunked_keys = [cmap_data.keys()[i::nthreads] for i in xrange(nthreads)]
    cmap_file_path_ext = os.path.splitext(cmap_file)
    for i in range(min(len(chunked_keys),nthreads)):
        curr_fname = cmap_file_path_ext[0] + "_chunk_" + str(i) + cmap_file_path_ext[1]
        with open(curr_fname,'w') as outfile:
            for l in cmap_header:
                outfile.write(l)

            for k in chunked_keys[i]:
                for l in cmap_data[k]:
                    outfile.write(l)

        cmap_chunked_fnames.append(curr_fname)

    return cmap_chunked_fnames

#Call SegAligner
class SegAlignerThread(threading.Thread):
    def __init__(self, threadID, contigs_file, ref_file, arg_list):
        threading.Thread.__init__(self)
        self.threadID = str(threadID)
        self.contigs_file = contigs_file
        self.ref_file = ref_file
        self.arg_list = arg_list

    def run(self):
        print "Starting " + self.name
        argstring = " ".join(self.arg_list)
        subprocess.call(" ".join(["SegAligner/SegAligner",self.ref_file,self.contigs_file,argstring]), shell=True)
        print "Finished thread " + self.threadID
        return 1

def read_scoring_file(scoring_file,scoring_dict):
    with open(scoring_file) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.rstrip().rsplit("\t")
            fields[0] = int(fields[0])
            fields[2] = float(fields[2])           
            scoring_dict[fields[0]].append((fields[1],fields[2]))

def get_relevant_contigs(seg_cutoffs,scoring_dict):
    rel_contigs = set()
    for seg_id,vals in scoring_dict.iteritems():
        S_cutoff = seg_cutoffs[seg_id]
        sig_hits_ind = bisect.bisect([x[1] for x in vals],S_cutoff)
        sig_vals = vals[sig_hits_ind:]
        for i in sig_vals:
            rel_contigs.add(i[0])

    return rel_contigs

def compute_e_value(scoring_dict,ref_file,p_val):
    #estimates an E-value (see Altschul and Karlin)
    ref_cmaps = parse_cmap(ref_file)
    for i in scoring_dict:
        scoring_dict[i] = sorted(scoring_dict[i],key=lambda x: x[1])

    m = reduce(lambda x,y:x+len(contig_cmaps[y]), contig_cmaps,0)

    seg_cutoffs = {}

    for i in scoring_dict:
        n = len(ref_cmaps[str(abs(i))])
        all_score_vals = [x[1] for x in scoring_dict[i]]
        #prune the outlier hits using Tukey's fences
        quarts = np.percentile(all_score_vals,[25,50,75])
        left_side = quarts[0] - 1.5*(quarts[2] - quarts[0])
        right_side = quarts[2] + 1.5*(quarts[2] - quarts[0])

        #get rid of bottom 50% of data
        left_ind = len(all_score_vals)/2 + 1

        #TUKEY
        # right_ind = bisect.bisect_left(all_score_vals,right_side)

        #top50
        right_ind = -50

        #outlier-removed data
        score_vals = all_score_vals[left_ind:right_ind]        
        e_vals = np.log(range(1,len(score_vals)+1))

        #linear regression
        slope,intercept,_,_,_ = linregress(score_vals[::-1],e_vals)

        #solve for k
        K = np.exp(intercept - np.log(m*n))

        #compute E_value corresponding to p_value cutoff
        E_cutoff = -np.log(1-p_val)
        S_cutoff = -np.log(E_cutoff/(K*m*n))/(-slope)
        seg_cutoffs[i] = S_cutoff
        print i,"smallest E-value found: ",(K*m*n)*np.exp(slope*all_score_vals[-1])

    return seg_cutoffs

def get_unaligned_segs(aln_path,aln_flist):
    contig_aln_dict = defaultdict(list)
    contig_unaligned_regions = defaultdict(list)
    for f in aln_flist:
        #parse the alnfile
        contig_id,aln_struct = parse_seg_alignment_file(aln_path + f)
        contig_aln_dict[contig_id].append(aln_struct)

    #extract aligned regions
    for c_id,a_struct_list in contig_aln_dict.iteritems():
        contig_free_lab_set = set(range(len(contig_cmaps[c_id])))
        for a_struct in a_struct_list:
            curr_c_range_tup = a_struct[2]
            curr_lab_range_set = set(range(curr_c_range_tup[0],curr_c_range_tup[1]+1))
            contig_free_lab_set-=curr_lab_range_set

        #extract unaligned regions using aligned region labels
        sorted_contig_free_labs = sorted(list(contig_free_lab_set))
        chunk_list = []
        prev = -2
        for i in sorted_contig_free_labs:
            if i - prev > 1:
                chunk_list.append([])
            
            chunk_list[-1].append(i)
            prev = i

        #check if unaligned region is large
        for l in chunk_list:
            unaligned_size = contig_cmaps[c_id][l[-1]+1] - contig_cmaps[c_id][l[0]+1]
            if l[-1] - l[0] > unaligned_label_cutoff and unaligned_size > unaligned_size_cutoff:
                #add the unaligned region
                contig_unaligned_regions[c_id].append(l)


    return contig_unaligned_regions

#write the unaligned regions, do some bookkeeping, mapping contig_id to some set of labels.
#return a mapping of the new contig_id and labels involved
def write_unaligned_cmaps(contig_unaligned_regions,output_prefix,enzyme):
    max_contig_id = max([int(x) for x in contig_cmaps.keys()])
    unaligned_region_contig_dict = {}
    unaligned_contig_id_dict = {}
    nmaps = sum([len(x) for x in contig_unaligned_regions.values()])
    cmap_header = "# CMAP File Version: 0.1\n# Label Channels:  1\n# Nickase Recognition Site 1:    " + enzyme + "\n# Enzyme1:  " + enzyme + "\n# Number of Consensus Nanomaps: " + str(nmaps) + "\n#h CMapId ContigLength    NumSites    SiteID  LabelChannel    Position    StdDev  Coverage    Occurrence\n#f int  float   int int int float   float   int int\n"
    unaligned_region_filename = output_prefix + "contig_unaligned_regions.cmap"
    with open(unaligned_region_filename,'w') as outfile:
        outfile.write(cmap_header)
        total_unaligned_segs = 0
        for c_id,all_regions in contig_unaligned_regions.iteritems():
            for l in all_regions:
                total_unaligned_segs+=1
                try:
                    map_end_pos = contig_cmaps[c_id][l[-1]+2]
                except KeyError:
                    map_end_pos = contig_lens[c_id]

                pos_l = [contig_cmaps[c_id][x+1] for x in l]
                #must be int
                # curr_id = c_id + "_" + str(l[0]+1) + "_" + str(l[-1]+1)
                curr_id = str(max_contig_id+total_unaligned_segs)
                unaligned_region_contig_dict[curr_id] = {str(ind+1):str(x+1) for ind,x in enumerate(l)}
                unaligned_contig_id_dict[curr_id] = c_id

                p0 = pos_l[0]
                map_len_str = str(map_end_pos - p0)
                for ind,i in enumerate(pos_l):
                    out_list = [curr_id,map_len_str,str(len(l)),str(ind+1),"1",str(i-p0),"1.0","1","1"]
                    outfile.write("\t".join(out_list) + "\n")

                out_list = [curr_id,map_len_str,str(len(l)),str(ind+1),"0",map_len_str,"1.0","1","1"]
                outfile.write("\t".join(out_list) + "\n")

    return unaligned_region_contig_dict,unaligned_region_filename,unaligned_contig_id_dict

#deprecated an unused function
def detections_to_graph(outfile,bpg_list):
    print "INSIDE", bpg_list
    for i in bpg_list:
        print i
        seg_size = int(i[2]) - int(i[1])
        outfile.write("sequence\t")
        outfile.write(i[0] + ":" + i[1] + "-\t")
        outfile.write(i[0] + ":" + i[2] + "+\t")
        outfile.write("2.0\t0\t" + str(seg_size) + "\t0\n")

def detections_to_key(outfile,keyfile_info):
    for i in keyfile_info:
        outfile.write(i[0] + "\t" + i[1] + ":" + i[2] + "-|" + i[1] + ":" + i[3] + "+\t0\n")

def rewrite_graph_and_CMAP(segs_fname,graphfile,bpg_list,enzyme):
    #read graph
    graphfile_lines = []
    with open(graphfile) as infile:
        for line in infile:
            graphfile_lines.append(line)

    new_graphfile = os.path.splitext(graphfile)[0] + "_includes_detected.txt"
    with open(new_graphfile, 'w') as outfile:
        for line in graphfile_lines:
            if line.startswith("BreakpointEdge:"):
                print "HERE"
                detections_to_graph(outfile,bpg_list)

            outfile.write(line)

    print "Creating new CMAP"
    subprocess.call("python graph_to_cmap.py -i " + new_graphfile + " -r ref_genomes/hg19.fa -e " + enzyme,shell=True)


def detections_to_seg_alignments(w_dir,aln_files,ref_file,unaligned_cid_d,unaligned_label_trans,id_start):
    #must indicate that the reference genome used (field in the head)
    ref_genome_cmaps = parse_cmap(ref_file)
    ref_genome_key_file = os.path.splitext(ref_file)[0] + "_key.txt"
    ref_genome_key_dict = parse_keyfile(ref_genome_key_file)

    seg_dir_count = defaultdict(int)
    bpg_list = []
    aln_num = 0
    print unaligned_cid_d.keys()
    print "aln_files len",len(aln_files)
    for f in aln_files:
        f_fields = os.path.splitext(f)
        f_fields = os.path.splitext(f)
        # outname = w_dir + f_fields[0] + "_corrected" + f_fields[1]
        head_list = []
        aln_field_list = []
        with open(w_dir + f) as infile:
            for i in range(3):
                head_list.append(infile.next())

            for line in infile:
                aln_field_list.append(line.rstrip().rsplit("\t"))

        #get direction 
        meta_list = []
        contig_id = head_list[1].rsplit("\t")[0][1:-1]
        contig_dir = head_list[1].rsplit("\t")[0][-1]
        try:
            trans_contig_id = unaligned_cid_d[contig_id]
        except KeyError:
            continue

        aln_num+=1
        #set up the breakpoint graph stuff
        ref_genome_id = aln_field_list[0][0]
        chromID = ref_genome_key_dict[ref_genome_id][0]
        p1 = int(ref_genome_cmaps[ref_genome_id][int(aln_field_list[0][2])])
        p2 = int(ref_genome_cmaps[ref_genome_id][int(aln_field_list[-1][2])])
        if p1 > p2:
            temp = p2
            p2 = p1
            p1 = temp

        #padding with 10 to see if it helps
        bpg_list.append((chromID,str(p1-10),str(p2+10)))

        # meta_list.append(trans_contig_id)
        unique_id = str(id_start + aln_num)
        meta_list.append(unique_id)
        if contig_dir == "-":
            aln_field_list = aln_field_list[::-1]

        meta_list[0]+=contig_dir
        meta_list.append(head_list[1].rsplit("\t")[1])
        meta_list.append("False")
        meta_string = "#" + "\t".join(meta_list) + "\n"

        # reverse the reference sequence if "-"
        isRev = "_r" if contig_dir == "-" else ""

        outname = "segalign_" + trans_contig_id + "_" + unique_id + "_rg" + isRev + "_aln.txt"
        with open(w_dir + outname,'w') as outfile:
            outfile.write(head_list[0])
            outfile.write(meta_string)
            outfile.write(head_list[2])
            curr_score = 0.0
            lchange = 0
            first = min(int(aln_field_list[0][2]),int(aln_field_list[-1][2]))
            for l in aln_field_list:
                trans_lab = unaligned_label_trans[contig_id][l[3]]
                new_l = [trans_contig_id,unique_id,trans_lab,str(int(l[2])-first+1)] + l[4:7] + [str(curr_score),str(lchange)]
                curr_score+=float(l[-1])
                lchange = l[-1]
                outfile.write("\t".join(new_l) + "\n")


        seg_dir_count[(trans_contig_id,contig_dir)]+=1

    return bpg_list

def run_SegAligner(contig_flist,segs,arg_list,nthreads):
    threadL = []
    for i in range(min(nthreads,len(contig_flist))):
        threadL.append(SegAlignerThread(i,contig_flist[i],segs,arg_list))
        threadL[i].start()

    for t in threadL:
        t.join()

#main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Banded DP alignment with scoring optimizations for BioNano data")
    parser.add_argument("-s", "--segs", help="reference genome cmap file. Key file must be present in same directory and be named [CMAP]_key.txt", required=True)
    parser.add_argument("-c", "--contigs", help="contigs cmap file", required=True)
    parser.add_argument("-g", help="AA breakpoint graph file",required=True)
    parser.add_argument("-r", "--ref", help="reference genome cmap file")
    parser.add_argument("-p", "--p_value", help="p-value threshold for alignments", type=float, default=0.05)
    parser.add_argument("-o", "--output_prefix", help="output filename prefix (assumes working directory & ref name unless otherwise specified")
    parser.add_argument("-t", help="number of threads to use (default 4)", type=int, default=4)
    parser.add_argument("-e", help="labeling enzyme", default="BspQI")
    args = parser.parse_args()


    if not args.output_prefix:
        args.output_prefix = os.path.splitext(args.segs)[0]
        print "results in " + args.output_prefix

    if not args.output_prefix.endswith("/"): args.output_prefix+="/"

    a_dir = args.output_prefix + "alignments/"
    if not os.path.exists(a_dir):
        os.makedirs(a_dir)


    min_map_len = 4
    nthreads = args.t

    contig_cmaps = parse_cmap(args.contigs)
    seg_cmaps = parse_cmap(args.segs)
    contig_lens = get_cmap_lens(args.contigs)
    
    #get chunked cmap
    chunked_cmap_files = chunk_cmap(args.contigs,nthreads)

    print "\nGetting distribution of best alignment scores for reference with contigs"
    print "Starting scoring at " + time.ctime(time.time())
    #CONTIG SCORING
    run_SegAligner(chunked_cmap_files,args.segs,[str(min_map_len),"-scoring","-prefix=" + args.output_prefix],nthreads)
    print "Scoring finished at " + time.ctime(time.time()) + "\n"

    #get the scoring thresholds.
    scoring_dict = defaultdict(list)
    flist = os.listdir(args.output_prefix)
    for f in flist:
        if f.startswith("scores"):
            read_scoring_file(args.output_prefix + f,scoring_dict)

    seg_cutoffs = compute_e_value(scoring_dict,args.segs,args.p_value)
    rel_contigs = get_relevant_contigs(seg_cutoffs,scoring_dict)
    #write the cutoffs and the contig list to files
    with open(args.output_prefix + "scoring_thresholds.txt",'w') as outfile:
        for key,val in seg_cutoffs.iteritems():
            outfile.write(str(key) + "\t" + str(val) + "\n")

    with open(args.output_prefix + "relevant_contigs.txt",'w') as outfile:
        for i in rel_contigs:
            outfile.write(str(i) + "\n")

    print "Doing alignments "
    arg_list = [str(min_map_len),"-prefix=" + a_dir,
    "-score_thresh=" + args.output_prefix + "scoring_thresholds.txt", 
    "-contig_list=" + args.output_prefix + "relevant_contigs.txt"]
    # #CONTIG SEG ALIGNMENTS
    run_SegAligner(chunked_cmap_files,args.segs,arg_list,nthreads)
    print "Alignments finished at " + time.ctime(time.time()) + "\n"

    #extract unaligned regions
    aln_flist = [x for x in os.listdir(a_dir) if x.startswith("segalign") and "flipped" not in x and "rg" not in x]
    #extract the unaligned regions.
    contig_unaligned_regions = get_unaligned_segs(a_dir,aln_flist)
    unaligned_label_trans,unaligned_region_filename,unaligned_cid_d = write_unaligned_cmaps(contig_unaligned_regions,args.output_prefix,args.e)

    #SCORING OF UNALIGNED REGIONS
    chunked_ref_files = []
    if contig_unaligned_regions:
        print "doing unaligned region detection"
        chunked_ref_files = chunk_cmap("ref_genomes/hg19_" + args.e + ".cmap",nthreads,500)
        run_SegAligner(chunked_ref_files,unaligned_region_filename,["-detection","-prefix=" + args.output_prefix],nthreads)
        
        #read the detected scores files
        #get the scoring thresholds.
        scoring_dict = defaultdict(list)
        flist = os.listdir(args.output_prefix)
        for f in flist:
            if f.startswith("detection_scores"):
                read_scoring_file(args.output_prefix + f,scoring_dict)

        #compute an e-value
        seg_cutoffs = compute_e_value(scoring_dict,unaligned_region_filename,args.p_value)
        rel_contigs = get_relevant_contigs(seg_cutoffs,scoring_dict)
        #write the cutoffs and the contig list to files
        with open(args.output_prefix + "detection_scoring_thresholds.txt",'w') as outfile:
            for key,val in seg_cutoffs.iteritems():
                outfile.write(str(key) + "\t" + str(val) + "\n")

        with open(args.output_prefix + "detection_relevant_contigs.txt",'w') as outfile:
            for i in rel_contigs:
                outfile.write(str(i) + "\n")

        #call alignment
        print "Doing alignments "
        arg_list = [str(min_map_len),"-prefix=" + a_dir,
        "-score_thresh=" + args.output_prefix + "detection_scoring_thresholds.txt", 
        "-contig_list=" + args.output_prefix + "detection_relevant_contigs.txt", "-local","-alnref"]
        #CONTIG UNALIGNED REGION ALIGNMENTS
        run_SegAligner(chunked_ref_files,unaligned_region_filename,arg_list,nthreads)
        print "Alignments finished at " + time.ctime(time.time()) + "\n"
        aln_flist = os.listdir(a_dir)
        with open(args.g) as infile:
            index_start = sum(1 for _ in infile if _.startswith("sequence"))

        bpg_list = detections_to_seg_alignments(a_dir,aln_flist,"ref_genomes/hg19_" + args.e + ".cmap",unaligned_cid_d,unaligned_label_trans,index_start)
        print "Found new segments, re-writing graph and CMAP"
        rewrite_graph_and_CMAP(args.segs,args.g,bpg_list,args.e)

    else:
        print "No large unaligned regions found on segment-aligned contigs"


    #remove the chunked cmaps
    print "removing temporary files"
    # subprocess.call("rm " + a_dir + "*flipped_aln.txt", shell=True)
    for i in chunked_cmap_files + chunked_ref_files:
        subprocess.call("rm " + i, shell=True)

    print "Completed " + time.ctime(time.time()) + "\n"










