from collections import defaultdict

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

#TODO: REFACTOR
#connect overlapping contig graphs
def find_intercontig_edges(scaffold_paths,contig_graphs):
    #take the scaffold heaviest paths and convert unique id to named segments for
    #overlap identification
    ##put all the graphs together, no intercontig edges present yet
    scaffold_paths_named = defaultdict(list)
    for c_id,path_list in scaffold_paths.iteritems():
        G = contig_graphs[c_id]
        for hp_ids in path_list:
            curr_path = [(G.node_id_lookup[i].seg_id,G.node_id_lookup[i].direction) for i in hp_ids]
            scaffold_paths_named[c_id].append(curr_path)

    suffix_d = {}
    prefix_d_fwd = {}
    prefix_d_rev = {}
    suffix_d_rev = {}
    for c_id, path_list in scaffold_paths_named.iteritems():
        for path_ind,curr_path in enumerate(path_list):
            path_len = len(curr_path)
            #Construct lookup for suffixes and map suffixes to graph_id
            for ind in range(path_len):
                curr_suffix = ",".join([str(-1*int(x[0])) if x[1] == "-" else x[0] for x in curr_path[ind:]])
                if curr_suffix not in suffix_d:
                    suffix_d[curr_suffix] = []

                suffix_d[curr_suffix].append((c_id,ind,path_ind))

            #prefixes in the forward direction
            for ind in range(1,path_len):
                curr_prefix = ",".join([str(-1*int(x[0])) if x[1] == "-" else x[0] for x in curr_path[:ind]])
                if curr_prefix not in prefix_d_fwd:
                    prefix_d_fwd[curr_prefix] = []

                prefix_d_fwd[curr_prefix].append((c_id,ind-1,path_ind)) #up to but not including ind

            #prefixes in the reverse direction
            for ind in range(1,path_len):
                curr_num_end_seq = [-1*int(x[0]) if x[1] == "-" else int(x[0]) for x in curr_path[ind:]]
                rev_num_seq = [str(-1*x) for x in curr_num_end_seq[::-1]]
                rev_prefix = ",".join(rev_num_seq)
                if not rev_prefix in prefix_d_rev:
                    prefix_d_rev[rev_prefix] = []

                prefix_d_rev[rev_prefix].append((c_id,ind,path_ind))

            #suffixes in the reverse direction
            for ind in range(1,path_len):
                curr_num_end_seq = [-1*int(x[0]) if x[1] == "-" else int(x[0]) for x in curr_path[:ind]]
                rev_num_seq = [str(-1*x) for x in curr_num_end_seq[::-1]]
                rev_prefix = ",".join(rev_num_seq)
                if not rev_prefix in suffix_d_rev:
                    suffix_d_rev[rev_prefix] = []

                suffix_d_rev[rev_prefix].append((c_id,ind-1,path_ind))


    s_p_matches = []
    suffix_set = set(suffix_d.keys())
    print suffix_set
    #intersect fwd set with suffix
    fwd_prefix_set = set(prefix_d_fwd.keys())
    fwd_intersect_keys = suffix_set.intersection(fwd_prefix_set)
    print fwd_prefix_set

    #intersect rev set with suffix
    rev_prefix_set = set(prefix_d_rev.keys())
    rev_intersect_keys = suffix_set.intersection(rev_prefix_set)
    print rev_prefix_set
    #

    rev_suffix_set = set(suffix_d_rev.keys())
    fwd_rev_suff_intersect_keys = rev_suffix_set.intersection(fwd_prefix_set)

    #make a set of new edges
    intercontig_edges = set()

    #connect matches with an inter-contig edge.
    #the tuple
    pairings = [(fwd_intersect_keys,prefix_d_fwd,0), (rev_intersect_keys,prefix_d_rev,-1), ()]

    #
    added_edges = set()
    for intersect_keys,prefix_d_curr,end_ind in zip([fwd_intersect_keys,rev_intersect_keys],[prefix_d_fwd,prefix_d_rev],[0,-1]):
        for i in intersect_keys:
            suffix_tup_list,prefix_tup_list = suffix_d[i],prefix_d_curr[i]
            for prefix_tup in prefix_tup_list:
                for suffix_tup in suffix_tup_list:
                    suffix_cid,suffix_ind,s_path_ind = suffix_tup
                    prefix_cid,prefix_ind,p_path_ind = prefix_tup
                    if prefix_cid == suffix_cid and prefix_d_curr == prefix_d_rev:
                        continue

                    Gs = contig_graphs[suffix_cid]
                    Gp = contig_graphs[prefix_cid]
                    hp_ids_suff = scaffold_paths[suffix_cid][s_path_ind]
                    hp_ids_pre = scaffold_paths[prefix_cid][p_path_ind]
                    s = Gs.node_id_lookup[hp_ids_suff[suffix_ind]]
                    t = Gp.node_id_lookup[hp_ids_pre[end_ind]]
                    if (s.n_id,t.n_id) in added_edges:
                        continue

                    new_edge = segment_edge(s,t,False)
                    new_edge.intercontig = True
                    if prefix_d_curr == prefix_d_rev:
                        new_edge.orientation_flip = True

                    print new_edge.edge_to_string(),end_ind
                    intercontig_edges.add(new_edge)
                    added_edges.add((s.n_id,t.n_id))
                    added_edges.add((t.n_id,s.n_id))

    for i in fwd_rev_suff_intersect_keys:
        suffix_tup_list,prefix_tup_list = suffix_d_rev[i],prefix_d_fwd[i]
        for prefix_tup in prefix_tup_list:
            for suffix_tup in suffix_tup_list:
                suffix_cid,suffix_ind,s_path_ind = suffix_tup
                prefix_cid,prefix_ind,p_path_ind = prefix_tup
                if prefix_cid == suffix_cid:
                    continue

                Gs = contig_graphs[suffix_cid]
                Gp = contig_graphs[prefix_cid]
                hp_ids_suff = scaffold_paths[suffix_cid][s_path_ind]
                hp_ids_pre = scaffold_paths[prefix_cid][p_path_ind]
                a = Gp.node_id_lookup[hp_ids_pre[0]]
                b = Gs.node_id_lookup[hp_ids_suff[suffix_ind]]
                if (a.n_id,b.n_id) in added_edges:
                    continue

                new_edge = segment_edge(a,b,False)
                new_edge.intercontig = True
                new_edge.orientation_flip = True

                print new_edge.edge_to_string(),end_ind
                intercontig_edges.add(new_edge)
                added_edges.add((a.n_id,b.n_id))
                added_edges.add((b.n_id,a.n_id))

    # for i in intercontig_edges:
    #     print(i.edge_to_string())

    return intercontig_edges