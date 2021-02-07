"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import sys
from collections import defaultdict

#scaffold alignment object
class SA_Obj(object):
    def __init__(self, contig_id, raw_aln_list):
        self.contig_id = contig_id
        self.seg_id = raw_aln_list[0]
        self.seg_endpoints = raw_aln_list[1]
        self.contig_endpoints = raw_aln_list[2]
        self.alignment_dir = raw_aln_list[3]
        self.aln_score = raw_aln_list[4]
        self.alignment = raw_aln_list[5]
        self.is_tip_aln = raw_aln_list[6]
        self.aln_id = "_".join([str(self.contig_id),str(self.seg_id)+self.alignment_dir,str(self.contig_endpoints[0]),str(self.contig_endpoints[1])])
        self.imputed_alignment = []
        self.is_detection_aln = False

    def aln_summary_to_string(self):
        return "seg_id: " + self.seg_id + " seg_labs: " + str(self.seg_endpoints) + " contig_labs: " + str(self.contig_endpoints) + " dir: " + self.alignment_dir + " aln_score: " + str(self.aln_score)

    def aln_to_string_list(self):
        return [str(x) for x in [self.seg_id,self.seg_endpoints,self.contig_endpoints,self.alignment_dir]]


#node class for alignment
class segment_node(object):
    def __init__(self,contig_id,aln_obj,imputed=False):
        self.aa_e = None
        self.aln_obj = aln_obj
        self.contig_id = contig_id
        self.seg_id = self.aln_obj.seg_id
        self.n_id = self.aln_obj.aln_id
        self.direction = self.aln_obj.alignment_dir
        self.imputed = imputed

    def node_to_string(self):
        return str((self.seg_id,self.contig_id,self.direction,self.n_id))

#node class for connections between alignments
class segment_edge(object):
    def __init__(self,s,t,f_status,imputed=False):
        self.s = s
        self.t = t
        self.forbidden = f_status
        self.junction_score = 0
        self.imputed = imputed
        self.intercontig = False
        self.gap = False
        self.heaviest_path_edge = False
        self.orientation_flip = False
        self.suboptimal = False

    def edge_to_string(self):
        return str((self.s.seg_id + self.s.direction,self.t.seg_id + self.t.direction,self.s.contig_id,self.t.contig_id,self.intercontig,self.imputed,self.orientation_flip))

#graph class for aligned nodes (the scaffold graph)
class contig_alignment_graph(object):
    def __init__(self):
        # self.ordered_node_list = ordered_node_list
        self.nodes = set()
        self.edges = set()
        self.adj_fwd_dict = defaultdict(list)
        self.adj_rev_dict = defaultdict(list)
        self.node_id_lookup = {}
        self.edge_lookup = {}
        self.weights = {}

    def construct_edge_lookup(self):
        self.edge_lookup = {(i.s.n_id,i.t.n_id):i for i in self.edges}

    def construct_adj_fwd_dict(self):
        for i in self.edges:
            if not i.forbidden:
                self.adj_fwd_dict[i.s.n_id].append(i.t.n_id)
    
    def construct_adj_rev_dict(self):
        for i in self.edges:
            if not i.forbidden:
                self.adj_rev_dict[i.t.n_id].append(i.s.n_id)

    #this is a weird way to solve it
    def construct_weights(self):
        for i in self.edges:
            if not i.forbidden:
                self.weights[(i.s.n_id,i.t.n_id)] = i.junction_score + i.t.aln_obj.aln_score
                self.weights[(i.t.n_id,i.s.n_id)] = i.junction_score + i.t.aln_obj.aln_score

    def construct_node_id_lookup(self):
        self.node_id_lookup = {i.n_id:i for i in self.nodes}

#make graph from alignments, does not consider overlapping contigs
def make_contig_aln_graph(aln_obj_list,contig_id,long_gap_length,allowed_overlap=1,cmap_id_to_edge={},contig_cmap={}):
    G = contig_alignment_graph()
    match_AA = True if cmap_id_to_edge else False
    #sort align list by startpoint
    raw_sorted_aln_l = sorted(aln_obj_list,key=lambda x: x.contig_endpoints[0])

    #make a list of sorted_aln_l nodes
    sorted_node_l = []
    sorted_aln_l = []
    for i in raw_sorted_aln_l:
        curr_node = segment_node(contig_id,i)
        if match_AA:
            try:
                curr_node.aa_e = cmap_id_to_edge[curr_node.seg_id]
            except KeyError:
                sys.stderr.write("Segment " + curr_node.seg_id + " not found in BPG\n")
                sys.stderr.write("Alignment files may not match to breakpoint graph.\n\n")

        sorted_node_l.append(curr_node)
        G.nodes.add(curr_node)
        sorted_aln_l.append(i)

    for ind_i, i in enumerate(sorted_aln_l[:-1]):
        lc_end = float('inf')
        curr_next = 1
        while ind_i + curr_next < len(sorted_aln_l):
            ind_j = ind_i + curr_next
            j = sorted_aln_l[ind_j]
            if j.contig_endpoints[0] > lc_end:
                break

            #forbidden
            if j.contig_endpoints[0] < i.contig_endpoints[1] - allowed_overlap: #ALLOWING OVERHANG 1 by default
                curr_edge = segment_edge(sorted_node_l[ind_i],sorted_node_l[ind_j],True)

            #not forbidden
            else:
                curr_edge = segment_edge(sorted_node_l[ind_i],sorted_node_l[ind_j],False)

                #if using om data
                if contig_cmap:
                    if contig_cmap[j.contig_endpoints[0]] - contig_cmap[i.contig_endpoints[1]] > long_gap_length:
                        curr_edge.gap = True

                #if using real coordinates
                else:
                    if j.contig_endpoints[0] - i.contig_endpoints[1] > long_gap_length:
                        curr_edge.gap = True

                if lc_end == float('inf'):
                    curr_is_tip_aln = sorted_node_l[ind_j].aln_obj.is_tip_aln
                    curr_is_detected = sorted_node_l[ind_j].aln_obj.is_detection_aln
                    if not curr_is_tip_aln and not curr_is_detected:
                        lc_end = j.contig_endpoints[1]

            G.edges.add(curr_edge)
            curr_next+=1

    G.ordered_node_list = sorted_node_l
    return G

#by default check if it is the alignment at the end of a contig
#right = True checks if alignment at the start of a contig
def is_end_aln(G,node_id,contig_cmap,left=False):
    curr = G.node_id_lookup[node_id]
    contig_endpoints = curr.aln_obj.contig_endpoints

    if left:
        dist_delta = contig_cmap[contig_endpoints[0]]
        lab_delta = contig_endpoints[0]

    else:
        dist_delta = contig_cmap[max(contig_cmap.keys())] - contig_cmap[contig_endpoints[1]]
        lab_delta = max(contig_cmap.keys()) - contig_endpoints[1]

    return (dist_delta < 125000 and lab_delta < 15)


def get_intercontig_edges(scaffold_paths,contig_graphs,contig_cmaps):
    prefix_f = defaultdict(list)
    prefix_r = defaultdict(list)
    suffix_f = defaultdict(list)
    suffix_r = defaultdict(list)
    scaffold_paths_named = defaultdict(list)
    prefix_suffix_list = []
    for c_id,path_list in scaffold_paths.iteritems():
        G = contig_graphs[c_id]
        for path_ind,hp_ids in enumerate(path_list):
            curr_path = [(G.node_id_lookup[i].seg_id,G.node_id_lookup[i].direction) for i in hp_ids]
            fwd_num_seq = ["-" + x[0] if x[1] == "-" else x[0] for x in curr_path]
            rev_num_seq = [str(-1*int(x)) for x in fwd_num_seq[::-1]]
            path_len = len(curr_path)

            #now get pres and suffs (if they're okay)
            #get all the prefixes, in both orientations
            if is_end_aln(G,hp_ids[0],contig_cmaps[c_id],left=True):
                for ind in range(1,path_len):
                    curr_prefix_f = ",".join(fwd_num_seq[:ind])
                    curr_prefix_r = ",".join(rev_num_seq[ind:])
                    prefix_f[curr_prefix_f].append((c_id,ind-1,path_ind))
                    prefix_r[curr_prefix_r].append((c_id,path_len - ind - 1,path_ind))

            #catch the case where a cycle is embedded in a single contig
            else:
                curr_prefix_f = fwd_num_seq[0]
                prefix_f[curr_prefix_f].append((c_id,0,path_ind))


            #get all the suffixes in similar fashion
            if is_end_aln(G,hp_ids[-1],contig_cmaps[c_id]):
                for ind in range(1,path_len):
                    curr_suffix_f = ",".join(fwd_num_seq[ind:])
                    curr_suffix_r = ",".join(rev_num_seq[:ind])
                    suffix_f[curr_suffix_f].append((c_id,ind,path_ind))
                    suffix_r[curr_suffix_r].append((c_id,path_len - ind,path_ind))

            else:
                curr_suffix_f = fwd_num_seq[-1]
                suffix_f[curr_suffix_f].append((c_id,path_len-1,path_ind))


    p_f_set = set(prefix_f.keys())
    p_r_set = set(prefix_r.keys())
    s_f_set = set(suffix_f.keys())
    s_r_set = set(suffix_r.keys())

    ##identify the overlaps
    #S+ -> P+
    s_p_intersect_keys = s_f_set.intersection(p_f_set)
    # print s_p_intersect_keys
    # print ""

    #S+ -> S-
    s_s_intersect_keys = s_f_set.intersection(s_r_set)
    # print s_s_intersect_keys
    # print ""

    #P- -> P+
    p_p_intersect_keys = p_r_set.intersection(p_f_set)
    # print p_p_intersect_keys
    # print ""

    #make a set of new edges
    intercontig_edges = set()
    added_edges = set()

    #iterate over the suffix prefix sets and check if there are valid overlaps
    param_list = [(s_p_intersect_keys,False,False,-1,suffix_f,prefix_f),(s_s_intersect_keys,True,True,-1,suffix_f,suffix_r),
    (p_p_intersect_keys,True,True,0,prefix_r,prefix_f)]
    for intersect_keys,disallow_self,orientation_flip,s_ind,source_d,dest_d in param_list:
        for i in intersect_keys:
            # print i
            source_list, dest_list = source_d[i],dest_d[i]
            pairings = [(s,t) for s in source_list for t in dest_list]
            # print pairings
            for source_tup,dest_tup in pairings:
                s_cid,_,s_path_ind = source_tup
                t_cid,t_ind,t_path_ind = dest_tup
                # print s_cid,t_cid
                if disallow_self and s_cid == t_cid:
                    # print "DS"
                    continue


                Gs = contig_graphs[s_cid]
                Gt = contig_graphs[t_cid]
                hp_ids_s = scaffold_paths[s_cid][s_path_ind]
                hp_ids_t = scaffold_paths[t_cid][t_path_ind]
                s_nid = hp_ids_s[s_ind]
                t_nid = hp_ids_t[t_ind]
                s = Gs.node_id_lookup[s_nid]
                t = Gt.node_id_lookup[t_nid]
                # print s_nid,t_nid

                #disallow connecting if it is an internal node and the destination is off-contig
                #also disallow if the two nodes are not appropriately ordered in thise case
                #this is to handle the interior cycle case
                s_is_end = is_end_aln(Gs,s_nid,contig_cmaps[s_cid]) or is_end_aln(Gs,s_nid,contig_cmaps[s_cid],left=True)
                t_is_end = is_end_aln(Gt,t_nid,contig_cmaps[t_cid]) or is_end_aln(Gt,t_nid,contig_cmaps[t_cid],left=True)
                if not (s_is_end and t_is_end):
                    if s_cid == t_cid:
                        if s.aln_obj.contig_endpoints[0] < t.aln_obj.contig_endpoints[-1]:
                            continue

                    elif not s_is_end and not t_is_end:
                            continue

                if (s.n_id,t.n_id) in added_edges:
                    continue

                new_edge = segment_edge(s,t,False)
                new_edge.intercontig = True
                if orientation_flip:
                    new_edge.orientation_flip = True

                intercontig_edges.add(new_edge)
                added_edges.add((s.n_id,t.n_id))
                added_edges.add((t.n_id,s.n_id))

    return intercontig_edges