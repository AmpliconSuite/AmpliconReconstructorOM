#!/usr/bin/env python

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import argparse
import copy
import datetime
import json
import math
import os
from subprocess import call
import sys

if sys.version_info < (3, 0):
    from Queue import Queue

else:
    from queue import Queue

from ContigAlignmentGraph import *
from bionanoUtil import *
from breakpoint_graph import *

# threshholds for linking contigs
long_gap_length = 400000  # long gap threshold between alignments
long_gap_cost = 0
max_search_depth = 64
max_impute_paths = 1024
max_paths_to_keep = 500
max_impute_edge_cc = 20


# uses a single list to keep the found paths to avoid having to flatten after returning
def dfs_path_find(t, curr_path, exp_length, curr_length, p_paths, c_count_d, last_edge):
    s = curr_path[-1]
    # check if path exceeds size difference constraints
    if not (curr_length - exp_length > min(10000 * (len(curr_path)), 25000)):
        if s.vid == t.vid:
            # check that is not too short
            if not (exp_length - curr_length > min(10000 * (len(curr_path)), 25000)):
                if len(curr_path) % 2 == 0:
                    p_paths.append(curr_path)

        if len(curr_path) < max_search_depth and len(p_paths) <= max_impute_paths:
            # get adj verts
            s_edges = s.elist
            for edge in s_edges:
                edge_rep = edge.__repr__()
                if last_edge == edge.edge_type or (
                        last_edge in ["concordant", "discordant"] and edge.edge_type in ["concordant", "discordant"]):
                    continue

                # must obey copy count
                if c_count_d[edge_rep] >= min(edge_cc[edge_rep], max_impute_edge_cc):
                    continue

                u = edge.neighbor(s)
                edge_len = 0
                if edge.edge_type == "sequence":
                    edge_len += (abs(s.pos - u.pos) + 1)

                # no direct loopbacks on short segs
                if len(curr_path) > 1:
                    if u.vid == curr_path[-2].vid and edge_len < 100:
                        continue

                c_count_d[edge_rep] += 1
                # recursive call
                dfs_path_find(t, curr_path + [u], exp_length, curr_length + edge_len, p_paths, c_count_d,
                              edge.edge_type)
                c_count_d[edge_rep] -= 1


# uses a single list to keep the found paths to avoid having to flatten after returning
def bfs_path_find(s, t, exp_length, seg_overhang_sum):
    p_paths = []
    bfs_terminate_count = max_search_depth * max_impute_paths
    curr_pop_count = 0
    path_queue = Queue()
    # Four things are stored in the each queue item
    # path list,  cc_dict,   path_len,   last_edge_type
    path_queue.put(([s], defaultdict(int), seg_overhang_sum, "sequence"))
    while not path_queue.empty() and curr_pop_count < bfs_terminate_count:
        curr_path, c_count_d, curr_length, last_edge = path_queue.get()
        curr_last_node = curr_path[-1]
        # check that it is not too long
        path_too_long = (curr_length - exp_length > min(10000 * (len(curr_path)), 25000))
        path_too_short = (exp_length - curr_length > min(10000 * (len(curr_path)), 25000))
        if not path_too_long and not path_too_short:
            if curr_last_node.vid == t.vid and len(curr_path) % 2 == 0:
                p_paths.append(curr_path)

        if not path_too_long and len(p_paths) <= max_impute_paths:
            s_edges = curr_last_node.elist
            for edge in s_edges:
                if last_edge == edge.edge_type or (
                        last_edge in ["concordant", "discordant"] and edge.edge_type in ["concordant", "discordant"]):
                    continue

                edge_rep = edge.__repr__()
                # must obey copy count
                if c_count_d[edge_rep] >= min(edge_cc[edge_rep], max_impute_edge_cc):
                    continue

                u = edge.neighbor(curr_last_node)
                edge_len = 0
                if edge.edge_type == "sequence":
                    edge_len += (abs(curr_last_node.pos - u.pos) + 1)

                # no direct loopbacks on short segs
                if len(curr_path) > 1:
                    if u.vid == curr_path[-2].vid and edge_len < 100:
                        continue

                c_count_d[edge_rep] += 1
                path_queue.put((curr_path + [u], copy.copy(c_count_d), curr_length + edge_len, edge.edge_type))

        elif len(p_paths) > max_impute_paths:
            logging.warning("BFS too large, stopped at " + str(curr_pop_count))
            break

        curr_pop_count += 1

    return p_paths


# method for checking if a better imputed path exists between two nodes on an edge
def path_alignment_correction(G, c_id, contig_cmap, impute=True):
    align_vertex_list = G.ordered_node_list
    es_to_remove = set()
    es_to_add = set()
    # for each edge, find all imputable paths between the two nodes
    for e in G.edges:
        if e.forbidden:
            continue

        i = e.s
        j = e.t

        if not e.s.aa_e or not e.t.aa_e:
            e.junction_score = 0
            continue

        if i.aln_obj.contig_endpoints[1] >= j.aln_obj.contig_endpoints[0]:
            logging.info("Tight junction at " + e.edge_to_string())
            e.junction_score = 1
            continue

        # do arithmetic on boundaries
        s_end = i.aln_obj.seg_endpoints[-1]
        t_start = j.aln_obj.seg_endpoints[0]
        s_last_aln_pos = vectorized_segs[i.seg_id][s_end - 1]
        t_first_aln_pos = vectorized_segs[j.seg_id][t_start - 1]

        # ensure order the vertices by position  (handle reverse)
        v1, v2 = i.aa_e.v1, i.aa_e.v2
        if v1.pos > v2.pos:
            v1 = i.aa_e.v2
            v2 = i.aa_e.v1

        # If reverse alignment, do switching of s and t directions
        s_id = i.seg_id
        if i.direction == "-":
            s_id += "_r"
            s_end = len(vectorized_segs[i.seg_id]) - s_end
            s = v1
            s_overhang = s_last_aln_pos

        else:
            s = v2
            s_overhang = vectorized_segs[i.seg_id][-1] - s_last_aln_pos

        # order the j vertices by position (handle reverse)
        v1, v2 = j.aa_e.v1, j.aa_e.v2
        if v1.pos > v2.pos:
            v1 = j.aa_e.v2
            v2 = j.aa_e.v1

        t_id = j.seg_id
        if j.direction == "-":
            t_id += "_r"
            t_start = len(vectorized_segs[j.seg_id]) - t_start
            t = v2
            t_overhang = vectorized_segs[j.seg_id][-1] - t_first_aln_pos

        else:
            t = v1
            t_overhang = t_first_aln_pos

        contig_start_label = i.aln_obj.contig_endpoints[-1]
        contig_end_label = j.aln_obj.contig_endpoints[0]
        # look up the expected distance
        contig_distance = contig_cmap[contig_end_label] - contig_cmap[contig_start_label]
        seg_overhang_sum = s_overhang + t_overhang

        possible_paths = []
        c_count_d = defaultdict(int)
        curr_path = [s]

        # don't look if there's a oversized gap or imputation is off
        if impute and not e.gap:
            logging.info("Searching for paths on edge " + s.__repr__() + " " + t.__repr__())
            dfs_path_find(t, curr_path, contig_distance, seg_overhang_sum, possible_paths, c_count_d, "sequence")

        elif e.gap:
            logging.info("[gap] Not searching for paths on edge " + s.__repr__() + " " + t.__repr__())

        if len(possible_paths) > max_impute_paths:
            logging.warning("Too many paths from DFS, attempting limited BFS")
            bfs_possible_paths = bfs_path_find(s, t, contig_distance, seg_overhang_sum)
            if len(bfs_possible_paths) <= max_impute_paths:
                possible_paths = bfs_possible_paths
                logging.info("Found " + str(len(possible_paths)) + " path(s) for " + s.__repr__() + " " + t.__repr__())
            else:
                possible_paths = []
                logging.warning("Too many BFS paths")

        else:
            logging.info("Found " + str(len(possible_paths)) + " path(s) for " + s.__repr__() + " " + t.__repr__())

        if [s, t] not in possible_paths:
            possible_paths.append([s, t])

        best_path = []
        best_score = float('-inf')

        for path in possible_paths:
            # get path score
            if e.gap:
                e.junction_score = long_gap_cost
                p_score = float('-inf')

            else:
                # make contig local cmap
                contig_start_pos = contig_cmap[contig_start_label]
                local_contig_cmap = []
                local_contig_rev_lookup = {}
                for ind, x in enumerate(range(contig_start_label, contig_end_label + 1)):
                    local_contig_cmap.append(contig_cmap[x] - contig_start_pos)
                    local_contig_rev_lookup[x] = ind + 1

                # get compound cmap containing the individual segments
                compound_cmap, compound_rev_lookup = path_to_cmaps(path, s_id, s_end, t_id, t_start)
                p_score, seg_aln_obj = path_score_from_SA_fitting_aln(compound_cmap, local_contig_cmap, c_id)
                e.junction_score = p_score

            if p_score > best_score:
                best_score = p_score
                best_path = path
                best_aln = seg_aln_obj
                best_compound_rev_lookup = compound_rev_lookup

            if path == [s, t]:
                s_t_score = p_score

        if len(best_path) > 2 and best_score > s_t_score:
            e.suboptimal = True
            es_to_add |= add_path(G, best_path, best_aln, i, j, best_score, contig_start_label, contig_end_label,
                                  c_id, best_compound_rev_lookup)

    G.edges |= es_to_add
    # for e in G.edges:
    #     print e.edge_to_string()

    logging.info("Finished imputation")


# makes CMAP from imputed path
# returns combined & partial cmaps
def path_to_cmaps(path, s_id, s_end, t_id, t_start):
    # 0-indexed label lookups
    compound_rev_lookup = defaultdict(dict)
    unlabeled_segs_starts = []
    s_start_pos = segs_cmaps[s_id][s_end]
    compound_cmap = [0.0]
    # add on s overhang
    for ind, i in enumerate(vectorized_segs[s_id][s_end:-1]):
        compound_cmap.append(i - s_start_pos)
        compound_rev_lookup[s_id][len(compound_cmap) - 1] = s_end + ind

    curr_endpoint = vectorized_segs[s_id][-1] - s_start_pos

    if len(path) > 2:
        lst = path[1:-1]
        segment_bounds = zip(*[lst[i::2] for i in range(2)])
        for i in segment_bounds:
            curr_id = aa_id_to_cmap_id(i)
            for ind, j in enumerate(vectorized_segs[curr_id][:-1]):
                compound_cmap.append(curr_endpoint + j)
                compound_rev_lookup[curr_id][len(compound_cmap) - 1] = ind + 1

            # handle case where segment has no labels
            if not vectorized_segs[curr_id][:-1]:
                unlabeled_segs_starts.append((curr_id, curr_endpoint))

            # update endpoint
            curr_endpoint += vectorized_segs[curr_id][-1]

    # add on t overhang
    for ind, i in enumerate(vectorized_segs[t_id][:t_start]):
        compound_cmap.append(curr_endpoint + i)
        compound_rev_lookup[t_id][len(compound_cmap) - 1] = ind + 1

    return compound_cmap, compound_rev_lookup


# map AA graph edge id to cmap ID
def aa_id_to_cmap_id(i):
    a_rep = i[0].__repr__()
    b_rep = i[1].__repr__()
    key = a_rep + "|" + b_rep
    try:
        curr_id = seg_to_cmap_id[key]
    except KeyError:
        sys.stdout.write("Error: unfound key in keyfile " + key + ". Wrong input files?\n")
        return None

    return curr_id


# match the cmap ID to its AA graph edge
def match_cmap_graph_edge(aa_graph):
    cmap_id_to_edge = {}
    for curr_edge in aa_graph.es.values():
        if curr_edge.edge_type == "sequence":
            seg_id = curr_edge.__repr__().replace("->", "|")
            try:
                cmap_id = seg_to_cmap_id[seg_id]
                cmap_id_to_edge[cmap_id] = curr_edge
            except KeyError:
                sys.stderr.write("ERROR: Couldn't locate segment " + seg_id + "\n")
                sys.exit()

    return cmap_id_to_edge


# change the imputed segment alignment label numbering back to normal
# returns an alignment formatted like bionanoUtil's method
def relabel_seg_aln(seg_aln_obj, seg_rev_lookup, contig_start_label, seg_id):
    # seg aln is in the tuple,score format
    relabeled_seg_aln = []
    seg_dir = "-" if seg_id.endswith("_r") else "+"
    for i in seg_aln_obj.alignment:
        curr_compound_seg_label = int(i["seg_label"]) - 1
        curr_contig_label = int(i["contig_label"]) + contig_start_label - 1
        if curr_compound_seg_label in seg_rev_lookup:
            lab_num = seg_rev_lookup[curr_compound_seg_label]
            if seg_dir == "-":
                lab_num = translate_reversed_label(vectorized_segs, seg_id, lab_num)

            aln_dict = dict(zip(["contig_label", "seg_label", "score", "seg_dir"],
                                [str(curr_contig_label), str(lab_num), i["score"], seg_dir]))
            relabeled_seg_aln.append(aln_dict)

    return relabeled_seg_aln


# places a new imputed path in the graph
def add_path(G, path, seg_aln_obj, i, j, p_score, c_start_l, c_end_l, contig_id, compound_rev_lookup):
    es_to_add = set()
    ordered_nodes = [i]
    # Make the entry i into a new graph node
    lst = path[1:-1]
    segment_bounds = zip(*[lst[x::2] for x in range(2)])

    # add the new path to the graph (nodes)
    for ind, x in enumerate(segment_bounds):
        curr_id = aa_id_to_cmap_id(x)
        # fix curr_id so not reverse named:
        curr_undirected_id = curr_id.split("_r")[0]
        aln_struct = []
        alignment_dir = '+' if x[0].pos < x[1].pos else '-'
        seg_endpoints = None
        contig_endpoints = (c_start_l, c_end_l)
        aln_score = -1

        new_aln_obj = SA_Obj(contig_id,
                             [curr_undirected_id, seg_endpoints, contig_endpoints, alignment_dir, aln_score, [], False])
        new_aln_obj.imputed_alignment = relabel_seg_aln(seg_aln_obj, compound_rev_lookup[curr_id], c_start_l, curr_id)
        new_node = segment_node(contig_id, new_aln_obj, imputed=True)
        new_node.n_id = new_node.n_id + ".imp." + str(ind)
        new_node.aa_e = cmap_id_to_edge[new_node.seg_id]  # references global variable cmap_id_to_edge
        # add the imputed alignment
        G.nodes.add(new_node)
        ordered_nodes.append(new_node)

    ordered_nodes.append(j)

    # add new edges to graph
    for x, y in zip(ordered_nodes[:-1], ordered_nodes[1:]):
        new_edge = segment_edge(x, y, False, imputed=True)
        es_to_add.add(new_edge)

    new_edge.junction_score = p_score

    return es_to_add


# takes compound_cmap: vector of posns, including first, last. This method must add on the length label
# sub_vect: the region of the contig to align with. This method must add on the length label
# c_id: contig id
# return a score, and aln_obj
def path_score_from_SA_fitting_aln(compound_cmap, sub_vect, c_id):
    # write the CMAPs to files
    compound_seg_fname = adir + "compound_segs.cmap"
    write_cmap_from_vector([compound_cmap + [(compound_cmap[-1] + 1)]], compound_seg_fname)

    sub_contig_fname = adir + "subcontig.cmap"
    write_cmap_from_vector([sub_vect + [(sub_vect[-1] + 1)]], sub_contig_fname)

    # call SA
    SA_SRC = os.environ["SA_SRC"]
    cmd_list = [SA_SRC + "/SegAligner", compound_seg_fname, sub_contig_fname, "-fitting", "-min_labs=0",
                "-prefix=" + adir + "SA_temp", "-gen=" + gen]
    with open(adir + "SA.stdout", 'w') as outfile:
        call(cmd_list, stdout=outfile)

    # read the files
    fitting_aln_file = adir + "SA_temp_1_1_fitting_aln.txt"
    a_c_id, a_list = parse_seg_alignment_file(fitting_aln_file)
    aln_obj = SA_Obj(a_c_id, a_list)
    aln_obj.contig_id = c_id
    call("rm " + fitting_aln_file, shell=True)
    return aln_obj.aln_score, aln_obj


# parse breakpoint graph file to get CN info
def get_edge_copy_counts(breakpoint_file):
    cc_dict = {}
    seq_edge_reps = set()
    with open(breakpoint_file) as infile:
        for line in infile:
            if line.rstrip():
                fields = line.rsplit()
                if line.startswith("sequence"):
                    curr_cc = max(math.ceil(float(fields[3])), 1.0)
                    e_rep = fields[1] + "->" + fields[2]
                    cc_dict[e_rep] = curr_cc
                    seq_edge_reps.add(e_rep)

                elif fields[0] in ["concordant", "discordant"]:
                    if fields[1] in cc_dict:
                        logging.error("Edge name collision. Using larger weight")
                        cc_dict[fields[1]] = max(cc_dict[fields[1]], math.ceil(float(fields[2])))

                    else:
                        cc_dict[fields[1]] = max(math.ceil(float(fields[2])), 1.0)

    return cc_dict, seq_edge_reps


def check_path_cc(G, path, cc_dict):
    path_edge_counts = defaultdict(int)
    last_node = G.node_id_lookup[path[-1][0]]
    p_seg_id = None
    p_contig_id = None
    for i in path:
        curr_node = G.node_id_lookup[i[0]]
        cn_repr = curr_node.aa_e.__repr__()
        if not (curr_node.seg_id == p_seg_id and curr_node.contig_id != p_contig_id):
            path_edge_counts[cn_repr] += 1

        p_seg_id = curr_node.seg_id
        p_contig_id = curr_node.contig_id

    # if path_is_circular(G,path):
    #     path_edge_counts[cn_repr]-=1
    if G.node_id_lookup[path[0][0]].seg_id == p_seg_id:
        path_edge_counts[cn_repr] -= 1

    # now find the path singletons with the max in the edge_cc, this is the scaling factor
    for i in range(1, max(path_edge_counts.values()) + 1):
        scale_val = i
        singletons = dict()
        for key, value in path_edge_counts.items():
            if value == scale_val:
                singletons[key] = cc_dict[key]

        if singletons:
            max_singleton = max(singletons.values()) / i
        else:
            continue

        # print i,max(singletons.values()),max_singleton

        # anything with > 1 path count, that is larger than the rounded scaled counts causes a fail.
        for seg_rep, value in path_edge_counts.items():
            # print seg_rep,value
            value -= 1
            scaled_cc = max(scale_val, round(cc_dict[seg_rep] / max_singleton))
            # print cc_dict[seg_rep]
            # print scaled_cc
            if scaled_cc < value and value >= i:
                return False

    return True


# graph to cytoscape js dict
def graphs_to_cytoscapejs_dict(G):
    graph = {}
    graph["nodes"] = []
    graph["edges"] = []
    for i in G.nodes:
        x = {"data": {}}
        p1, p2 = seg_key[i.seg_id][0].split("|")
        x_name = "Contig " + i.contig_id + ": Seg " + i.seg_id + " " + p1[:-1] + "->" + p2[:-1] + "  " + i.direction
        x["data"]["id"] = i.n_id
        x["data"]["name"] = x_name
        graph["nodes"].append(x.copy())

    for i in G.edges:
        x = {"data": {}}
        s_id = i.s.n_id
        t_id = i.t.n_id
        e_id = s_id + " -> " + t_id
        x["data"]["source"] = s_id
        x["data"]["target"] = t_id
        x["data"]["id"] = e_id
        x["data"]["lText"] = str(i.junction_score)

        if i.forbidden:
            x["data"]["lCol"] = "red"

        else:
            x["data"]["lCol"] = "grey"

        if i.imputed:
            x["data"]["lStyle"] = "dotted"
        else:
            x["data"]["lStyle"] = "solid"

        if i.gap:
            x["data"]["lStyle"] = "dashed"
            x["data"]["lCol"] = "orange"

        elif i.intercontig:
            x["data"]["lStyle"] = "solid"
            x["data"]["lCol"] = "green"

        if i.heaviest_path_edge:
            x["data"]["lCol"] = "blue"

        graph["edges"].append(x.copy())

    return graph


# put all nodes and edges into a single graph and call the scaffold linking
def construct_combined_graph(scaffold_graphs):
    G = contig_alignment_graph()
    for i in scaffold_graphs.values():
        # merge nodes and edges into G
        G.nodes |= i.nodes
        G.edges |= i.edges

    return G


# topological sort helper
def topological_sort_recursion(G, n_id, visited, topo_stack):
    visited[n_id] = True
    for i in G.adj_fwd_dict[n_id]:
        if not visited[i]:
            topological_sort_recursion(G, i, visited, topo_stack)

    topo_stack.append(n_id)


# topological sort of a graph
def topological_sort(G):
    topo_stack = []
    visited = {i: False for i in G.node_id_lookup}
    for i in G.node_id_lookup:
        if not visited[i]:
            topological_sort_recursion(G, i, visited, topo_stack)

    return topo_stack


# computes heaviest path in a DAG starting at s in linear time with a DP method
def get_scaffold_heaviest_path(contig_G, s, topo_sorted_ids, weight_dict):
    # init dists
    dist = {i: float("-inf") for i in topo_sorted_ids}
    dist[s.n_id] = s.aln_obj.aln_score
    backtracking_dict = {s.n_id: (None, 0.0)}

    while topo_sorted_ids:
        u_id = topo_sorted_ids.pop()
        if dist[u_id] != float("inf"):
            for i in contig_G.adj_fwd_dict[u_id]:
                new_score = weight_dict[(u_id, i)] + contig_G.edge_lookup[(u_id, i)].junction_score
                if dist[i] < dist[u_id] + weight_dict[(u_id, i)]:
                    dist[i] = dist[u_id] + weight_dict[(u_id, i)]
                    backtracking_dict[i] = (u_id, dist[i])

    # now backtrack to reconstruct path
    # get the largest value as the start
    curr = s.n_id
    max_weight = float("-inf")
    for i, value_tup in backtracking_dict.items():
        if value_tup[1] > max_weight:
            curr, max_weight = i, value_tup[1]

    heaviest_path = []
    while curr:
        heaviest_path.append(curr)
        curr = backtracking_dict[curr][0]

    return heaviest_path[::-1], max_weight


# recursive pathfinding for non-extendable paths
def path_recursion(G, u, visited, curr_path, paths, edge_dir, p_intercontig):
    if edge_dir == 1:
        unvisited_next = set([x for x in G.adj_fwd_dict[u]]) - visited
    else:
        unvisited_next = set([x for x in G.adj_rev_dict[u]]) - visited

    if not unvisited_next:
        # print "TN: " + path_to_string(G,curr_path,True)
        paths.append(curr_path)

    else:
        # print "PN: " + path_to_string(G,curr_path,True)
        # print "UN: " + str(unvisited_next)
        # print ""
        for v in unvisited_next:
            new_visited = copy.copy(visited)
            new_visited.add(v)
            n_tup = tuple([u, v][::edge_dir])
            curr_edge = G.edge_lookup[n_tup]
            if curr_edge.gap:
                paths.append(curr_path)
                continue

            elif curr_edge.suboptimal:
                continue

            new_edge_dir = edge_dir
            if curr_edge.orientation_flip:
                new_edge_dir *= -1

            if not (curr_edge.intercontig and p_intercontig):
                path_recursion(G, v, new_visited, curr_path + [(v, new_edge_dir)], paths, new_edge_dir,
                               curr_edge.intercontig)


# check if one path is entirely a subsequence of another
def check_LCS(path1, path2, downsample=False):
    x, y = len(path1), len(path2)
    M = [[0] * (y + 1) for i in range(x + 1)]
    for i in range(1, x + 1):
        for j in range(1, y + 1):
            if path1[i - 1] == path2[j - 1]:
                M[i][j] = M[i - 1][j - 1] + 1
            else:
                M[i][j] = max(M[i - 1][j], M[i][j - 1])

    if M[x][y] == min(x, y):
        return True

    # Handles long paths with minor modifications. Filtered later by score so better one is kept.
    elif downsample and min(x, y) > 12 and M[x][y] > min(x, y) - 3:
        return True

    return False


# Check if path violates copy count constraints from AA bgp
def filter_paths_by_cc(G, all_paths, edge_cc):
    cc_valid_paths = []
    for path in all_paths:
        if check_path_cc(G, path, edge_cc):
            cc_valid_paths.append(path)

    return cc_valid_paths


def get_segdir_seq(G, path):
    segdir_seq = [G.node_id_lookup[x[0]].seg_id + str(int(G.node_id_lookup[x[0]].direction + "1") * x[1]) for x in path]
    return segdir_seq


def check_rotations(G, kept, i, rev_i):
    i_segdir_seq = get_segdir_seq(G, i)
    rev_i_segdir_seq = get_segdir_seq(G, rev_i)
    i_circ = path_is_circular(G, i)
    for j in kept:
        j_circ = path_is_circular(G, j)
        j_segdir_seq = get_segdir_seq(G, j)
        if i_circ != j_circ:
            continue

        elif not i_circ:
            if len(j) > 1 or len(i) == 1:
                if check_LCS(i_segdir_seq, j_segdir_seq) or check_LCS(rev_i_segdir_seq, j_segdir_seq):
                    return True

        elif len(j) > 1:
            for rot_ind in range(len(j)):
                r_j = j[rot_ind:] + j[:rot_ind]
                r_j_segdir_seq = get_segdir_seq(G, r_j)
                if check_LCS(i_segdir_seq, r_j_segdir_seq) or check_LCS(rev_i_segdir_seq, r_j_segdir_seq):
                    return True

    return False


# return all paths from a list of paths which are not a sub-sequence
def filter_subsequence_paths(G, paths):
    kept = []
    contig_to_paths = defaultdict(list)
    logging.info("Sorting paths by weight")
    paths_sorted = sorted(paths, reverse=True, key=lambda x: get_path_weight(G, x))
    downsample = True if len(paths_sorted) > 35000 else False

    if downsample:
        logging.warning("Limiting search to top 35000 paths from " + str(len(paths_sorted)) + " original paths")
        paths_sorted = paths_sorted[:35000]

    for ind_i, i in enumerate(paths_sorted):
        if ind_i % 1000 == 1 and ind_i > 1:
            logging.info("Checked {}/{} paths, {} are still kept.".format(str(ind_i - 1), str(len(paths_sorted)),
                                                                          str(len(kept))))

        rev_i = [(x[0], -1 * x[1]) for x in i][::-1]
        if not check_rotations(G, kept, i, rev_i):
            kept.append(i)
            if downsample and len(kept) > max_paths_to_keep:
                break

    return kept


# get heaviest paths for each of the scaffold graphs
def all_unique_non_extendible_paths(G, edge_cc, scaffold_alt_paths, disable_CC_check=False):
    # construct all the intermediate nodes not to start at
    # (i.e. they are inside the heaviest path and not an endpoint)
    shp_interior_nodes = set()
    for c_id, path_list in scaffold_alt_paths.items():
        for path in path_list:
            for i in path[1:-1]:
                shp_interior_nodes.add(i)

    # iterate through nodes and recurse on the pseudo-directed graph to get the paths
    all_paths = []
    for i in [x.n_id for x in G.nodes if not x.imputed and x.n_id not in shp_interior_nodes]:
        paths = []
        path_recursion(G, i, {i}, [(i, 1)], paths, 1, True)
        all_paths.extend(paths)
        paths = []
        path_recursion(G, i, {i}, [(i, -1)], paths, -1, True)
        all_paths.extend(paths)

    # dump_paths_sorted = sorted(all_paths, reverse=True, key=lambda x: get_path_weight(G, x))
    # with open("dump.txt",'w') as outfile:
    #     for i in dump_paths_sorted:
    #         outfile.write(path_to_string(G,i,True) + "\n")

    logging.info("Total intial paths discovered: " + str(len(all_paths)))
    if not disable_CC_check:
        cc_paths = filter_paths_by_cc(G, all_paths, edge_cc)
        logging.info("Total CC filtered paths: " + str(len(cc_paths)))

    else:
        logging.info("SKIPPING CC CHECK, PER CLI ARGUMENT")
        cc_paths = all_paths

    # with open("ccdump.txt",'w') as outfile:
    #     for i in cc_paths:
    #         outfile.write(path_to_string(G,i,True) + "\n")

    ss_paths = filter_subsequence_paths(G, cc_paths)

    logging.info("Total final paths: " + str(len(ss_paths)))
    return ss_paths


# calculate the weight of a path
def get_path_weight(G, path):
    weight = G.node_id_lookup[path[0][0]].aln_obj.aln_score
    for s, t in zip(path[:-1], path[1:]):
        try:
            weight += G.weights[(s[0], t[0])]
        except KeyError:
            logging.error("BAD EDGE ON PATH " + path_to_string(G, path, True))

    if (path[-1][0], path[0][-1]) in G.edge_lookup:
        weight += G.edge_lookup[(path[-1][0], path[0][-1])].junction_score

    return weight


# use the directionality inferred by the path to set the final alignment directions
def get_final_direction(aln_dir, flipped):
    if flipped < 0:
        aln_dir = "+" if aln_dir == "-" else "-"

    return aln_dir


# return circularity and looping edge
def path_is_circular(G, path):
    seg_seq = [(G.node_id_lookup[i[0]], i[1]) for i in path]
    circular = False
    i = path[-1]
    e_tup = (path[0][0], path[-1][0]) if i[1] < 0 else (path[-1][0], path[0][0])
    if e_tup in G.edge_lookup:
        looping_edge = G.edge_lookup[e_tup]
        if not looping_edge.forbidden:
            circular = True

    return circular


# produce the AA-style cycle as a list
def path_to_cycle_list(G, path):
    # node_seq = [[G.node_id_lookup[i[0]],i[1]] for i in path]
    seg_seq = [(G.node_id_lookup[i[0]], i[1]) for i in path]
    circular = path_is_circular(G, path)

    # construct final cycle sequence and remove duplicate nodes caused by intercontig
    cycle_list = []
    contig_list = []
    aug_seg_seq = seg_seq + [seg_seq[0]]
    for ind, i in enumerate(aug_seg_seq[:-1]):
        oriented_segment_id = i[0].seg_id + get_final_direction(i[0].direction, i[1])
        if i[0].contig_id not in contig_list:
            contig_list.append(i[0].contig_id)

        # check if edge is intercontig, if it is then this thing is duplicated
        try:
            if i[1] < 0:
                curr_edge = G.edge_lookup[(aug_seg_seq[ind + 1][0].n_id, i[0].n_id)]

            else:
                curr_edge = G.edge_lookup[(i[0].n_id, aug_seg_seq[ind + 1][0].n_id)]

            if not curr_edge.intercontig:
                cycle_list.append(oriented_segment_id)

        except KeyError:
            if ind == len(seg_seq) - 1:
                cycle_list.append(oriented_segment_id)

            else:
                sys.stderr.write("Could not find expected interior edge in path at index " + str(ind))
                sys.stderr.write(path_to_string(G, path))

    if not circular:
        cycle_list = ["0+"] + cycle_list + ["0-"]

    elif circular and len(cycle_list) > 1 and cycle_list[-1] == cycle_list[0]:
        cycle_list.pop()

    return cycle_list, circular, contig_list


# write the path as an alignment file
def write_path_alignment(G, path, outname, weight):
    seg_seq, circular, _ = path_to_cycle_list(G, path)
    node_seq = [G.node_id_lookup[i[0]] for i in path]
    # print seg_seq,circular
    # print [x.n_id for x in node_seq]
    # print ""

    with open(outname, 'w') as outfile:
        outfile.write("#seg_seq\tmedian_aln_score\tmean_aln_score\ttotal_score\tcircular\n")
        outfile.write("#" + "\t".join([",".join(seg_seq), "0", "0", str(weight), str(circular)]) + "\n")
        outfile.write(
            "#contig_id\tseg_id\tcontig_label\tseg_label\tcontig_dir\tseg_dir\tseg_aln_number\tscore\tscore_delta\timputed\n")
        aln_ind = 0
        prev = None
        for ind, v in enumerate(node_seq):
            if not v.imputed:
                alignment, imputed = v.aln_obj.alignment, "0"
            else:
                alignment, imputed = v.aln_obj.imputed_alignment, "1"

            if ind != 0:
                lookup = (v.n_id, prev) if path[ind - 1][1] == -1 else (prev, v.n_id)
                if not G.edge_lookup[lookup].intercontig:
                    aln_ind += 1

            prev = v.n_id

            # #check if contig is negative direction:
            if path[ind][1] == 1:
                contig_dir = "+"
            else:
                contig_dir = "-"
                alignment = alignment[::-1]

            if ind == len(node_seq) - 1 and circular and node_seq[0].seg_id == node_seq[-1].seg_id:
                aln_ind = 0

            for i in alignment:
                seg_label = i["seg_label"]
                outlist = [v.contig_id, v.seg_id, i["contig_label"], seg_label, contig_dir, i["seg_dir"], str(aln_ind),
                           "0", "0", imputed]
                outfile.write("\t".join(outlist) + "\n")


# write the path as an AA cycles file
def write_path_cycles(G, paths, outname):
    with open(outname, 'w') as outfile:
        outfile.write("List of cycle segments\n")
        orig_segs = [x for x in segs_cmaps if "_r" not in x]
        for i in sorted(orig_segs, key=lambda x: int(x)):
            key_pos_string = seg_key[i][0][:-1]
            start, end = key_pos_string.split("-|")
            chrom, start = start.split(":")
            _, end = end.split(":")
            outfile.write("\t".join(["Segment", i, chrom, start, end]) + "\n")

        for ind, i in enumerate(paths):
            cycle_list, _, contig_list = path_to_cycle_list(G, i)
            outfile.write("Cycle=%d;Copy_count=1;Contigs=%s;Segments=%s\n" % (ind + 1, ",".join(contig_list),
                                                                              ",".join(cycle_list)))


def path_to_string(G, path, show_contig=False):
    line = ""
    for i in path:
        line += G.node_id_lookup[i[0]].seg_id
        line += G.node_id_lookup[i[0]].direction
        if show_contig:
            line += "(c_id:"
            line += G.node_id_lookup[i[0]].contig_id
            line += (")")

        line += ", "

    return line


def get_scaffold_heaviest_paths(contig_alignment_dict, impute, contig_cmaps):
    contig_graphs = {}
    scaffold_heaviest_paths = {}
    for ind, c_id in enumerate(contig_alignment_dict.keys()):
        logging.info("Path imputation, contig id: " + str(c_id))
        aln_obj_list = contig_alignment_dict[c_id]
        # put the contig segment alignments into the graph
        contig_graphs[c_id] = make_contig_aln_graph(aln_obj_list, c_id, long_gap_length,
                                                    cmap_id_to_edge=cmap_id_to_edge, contig_cmap=contig_cmaps[c_id])
        G_contig = contig_graphs[c_id]
        path_alignment_correction(G_contig, c_id, contig_cmaps[c_id], impute)
        G_contig.construct_edge_lookup()
        G_contig.construct_node_id_lookup()
        G_contig.construct_weights()
        G_contig.construct_adj_fwd_dict()
        best_path_weight = float("-inf")
        topo_sorted_ids = topological_sort(G_contig)
        for s in G_contig.nodes:
            heaviest_path, curr_weight = get_scaffold_heaviest_path(G_contig, s, copy.copy(topo_sorted_ids),
                                                                    G_contig.weights)
            if curr_weight > best_path_weight:
                best_path, best_path_weight = heaviest_path, curr_weight

            for s, t in zip(best_path[:-1], best_path[1:]):
                edge = G_contig.edge_lookup[(s, t)]
                edge.heaviest_path_edge = True

        scaffold_heaviest_paths[c_id] = (best_path, best_path_weight)

        # remove unused RG edges
        shp_node_set = set()
        for i in best_path:
            shp_node_set.add(i)

        unkept_rg_edges = set()
        for edge in G_contig.edges:
            if edge.s.aln_obj.is_detection_aln and edge.s.n_id not in shp_node_set:
                # print edge.s.n_id
                unkept_rg_edges.add(edge)

            elif edge.t.aln_obj.is_detection_aln and edge.t.n_id not in shp_node_set:
                # print edge.t.n_id
                unkept_rg_edges.add(edge)

        G_contig.edges -= unkept_rg_edges

        # DEBUGGING print
        logging.info("Heaviest path for contig " + c_id)
        logging.info(str([G_contig.node_id_lookup[i].seg_id for i in best_path]) + " " + str(best_path_weight) + "\n")

    return contig_graphs, scaffold_heaviest_paths


# add alternate "heaviest" paths, i.e. alternate paths formed by forbidden edge endpoints
def add_alternate_paths(contig_graphs, scaffold_heaviest_paths):
    # the set of paths considered for intercontig connections
    connectable_paths = {}
    for c_id, shp_tup in scaffold_heaviest_paths.items():
        if len(shp_tup[0]) == 1:
            connectable_paths[c_id] = [shp_tup[0], ]

        else:
            hp = shp_tup[0]
            first_node = hp[0]
            last_node = hp[-1]
            G_contig = contig_graphs[c_id]
            alt_first_nodes = [first_node, ]
            alt_last_nodes = [last_node, ]

            for i in [x for x in G_contig.nodes if not x.imputed]:
                for node, hp_i in zip([first_node, last_node], [hp[1], hp[-2]]):
                    if G_contig.node_id_lookup[hp_i].imputed:
                        continue

                    for node_pair in [(node, i.n_id), (i.n_id, node)]:
                        if node_pair in G_contig.edge_lookup:
                            curr_edge = G_contig.edge_lookup[node_pair]
                            if curr_edge.forbidden:
                                if (node, hp_i) in G_contig.edge_lookup and hp_i == hp[1]:
                                    if not G_contig.edge_lookup[(node, hp_i)].forbidden:
                                        alt_first_nodes.append(i.n_id)

                                elif (hp_i, node) in G_contig.edge_lookup and hp_i == hp[-2]:
                                    if not G_contig.edge_lookup[(hp_i, node)].forbidden:
                                        alt_last_nodes.append(i.n_id)

            alt_paths = []
            for s in alt_first_nodes:
                for t in alt_last_nodes:
                    alt_paths.append([s] + hp[1:-1] + [t])

            # print alt_paths
            connectable_paths[c_id] = alt_paths

    return connectable_paths


### MAIN ###
if __name__ == '__main__':
    # Parses the command line arguments
    parser = argparse.ArgumentParser(
        description="Corrects and extends alignment paths to produce OM contig/AA segment scaffolds")
    parser.add_argument("--adir", dest="adir", help="Directory where SegAligner alignments for the sample are located",
                        required=True)
    parser.add_argument("-c", "--contigs", help="contig cmap", required=True)
    parser.add_argument("-s", "--segs",
                        help="segments cmap file. Key file must be present in same directory and be named [CMAP]_key.txt",
                        required=True)
    parser.add_argument("-g", "--graph", help="AA graph file", required=True)
    parser.add_argument("--noImpute", help="Do not impute the paths", action='store_true', default=False)
    parser.add_argument("--outdir", type=str,
                        help="Destination for output files. Will create folder if it does not already exist",
                        required=True)
    parser.add_argument("--prefix", dest="samp_name", type=str, help="Filename prefix for output files.",
                        default="OMPF_out")
    parser.add_argument("--noConnect", action='store_true', help="Do not perform intercontig connection step")
    parser.add_argument("-i", "--instrument", choices=["Irys", "Saphyr"], required=True)
    parser.add_argument("--contig_subset", help="Use the following contig IDs only when reconstructing (default will "
                                                "use all contigs).", nargs='+', type=str, default=[])
    parser.add_argument("--disable_CC_check", help="Disable copy count ratio conforming check.", action='store_true',
                        default=False)

    args = parser.parse_args()

    impute = False if args.noImpute else True
    adir = args.adir
    if not adir.endswith("/"): adir += "/"

    # set up outdir
    outdir = args.outdir
    if not os.path.exists(outdir): os.makedirs(outdir)
    if not outdir.endswith("/"): outdir += "/"
    gen = "1" if args.instrument == "Irys" else "2"

    samp_name = args.samp_name
    outname = outdir + samp_name
    logging.basicConfig(filename=outname + "_run.log", level=logging.INFO)

    # read contig cmaps
    contig_cmaps_file = args.contigs
    contig_cmaps = parse_cmap(contig_cmaps_file, True)

    # read segs cmap
    segs_cmaps_file = args.segs
    segs_cmaps = parse_cmap(segs_cmaps_file, True)
    # print(segs_cmaps)

    # read segs key
    keyfile = os.path.splitext(args.segs)[0] + "_key.txt"
    if not os.path.exists(keyfile):
        sys.stderr.write(keyfile + " not found. Graph CMAP and *_key.txt must in same folder. Exiting.\n")
        sys.exit()
    seg_key = parse_keyfile(keyfile)

    # make reverse cmaps and reverse key
    add_full_reverse_cmaps(segs_cmaps, seg_key)
    seg_to_cmap_id = {v[0]: k for k, v in seg_key.items()}
    vectorized_segs = vectorize_cmaps(segs_cmaps)

    # reconstitute AA graph
    logging.info("Reconstituting breakpoint graph")
    breakpoint_file = args.graph
    breakpointG = breakpoint_graph(breakpoint_file)

    # get max copy count
    edge_cc, seq_edge_reps = get_edge_copy_counts(breakpoint_file)
    # unadj_edge_cc,seq_edge_reps = get_edge_copy_counts(breakpoint_file)
    # edge_cc = adjust_cc(unadj_edge_cc,seq_edge_reps)

    # match cmap to AA edges
    cmap_id_to_edge = match_cmap_graph_edge(breakpointG)
    # print(cmap_id_to_edge)
    logging.info("Matched ids to AA edges")

    # get aln files
    flist = os.listdir(adir)
    rel_files = []
    for i in flist:
        if i.endswith("_aln.txt") and not "_ref_" in i and not i.startswith("SA_temp"):
            rel_files.append(i)

    contig_alignment_dict = defaultdict(list)
    # parse each aln file
    logging.info("Parsing alignments")
    for i in rel_files:
        a_c_id, a_list = parse_seg_alignment_file(adir + i)
        seg_aln_obj = SA_Obj(a_c_id, a_list)
        c_id = seg_aln_obj.contig_id
        if len(args.contig_subset) > 0:
            if c_id not in args.contig_subset:
                continue

        seg_id = seg_aln_obj.seg_id
        if seg_id not in cmap_id_to_edge:
            sys.stderr.write(
                "Found segment id in alignments but not in breakpoint graph cmap. seg_id: " + str(seg_id) + "\n")
            sys.stderr.write(
                "This a non-fatal warning that the alignment files do not entirely match the breakpoint graph. Are you using the correct graph file?\n")
            continue

        if "_rg_" in i:
            seg_aln_obj.is_RG_aln = True

        contig_alignment_dict[c_id].append(seg_aln_obj)

    # Get the heaviest path per scaffold
    contig_graphs, scaffold_heaviest_paths = get_scaffold_heaviest_paths(contig_alignment_dict, impute, contig_cmaps)
    alt_paths = add_alternate_paths(contig_graphs, scaffold_heaviest_paths)

    # make intercontig edges
    G = construct_combined_graph(contig_graphs)

    if args.noConnect:
        logging.info("skipping contig connection step")
        intercontig_edges = set()
    else:
        intercontig_edges = get_intercontig_edges(alt_paths, contig_graphs, contig_cmaps)

    G.edges |= intercontig_edges
    G.construct_edge_lookup()
    G.construct_node_id_lookup()
    G.construct_weights()
    G.construct_adj_fwd_dict()
    G.construct_adj_rev_dict()

    # connect heaviest paths across contigs
    logging.info("Finding all non-extendible paths")
    all_paths = all_unique_non_extendible_paths(G, edge_cc, alt_paths, args.disable_CC_check)

    all_paths_weights = [get_path_weight(G, p) for p in all_paths]
    if args.noImpute: outname += "_noImpute"
    # write aligned path results
    # write discovered paths
    logging.info("Writing paths")
    if all_paths:
        sorted_all_paths, sorted_all_weights = zip(
            *sorted(zip(all_paths, all_paths_weights), key=lambda x: x[1], reverse=True))
        for ind, i in enumerate(sorted_all_paths):
            fname = "%s_path_%d_aln.txt" % (outname, ind + 1)
            write_path_alignment(G, i, fname, sorted_all_weights[ind])

        write_path_cycles(G, sorted_all_paths, outname + "_paths_cycles.txt")

    # write SHPs
    flattened_shps = [p[0] for p in scaffold_heaviest_paths.values()]
    flattened_directed_shps = []
    for p in flattened_shps:
        cl = []
        for n in p:
            cl.append((n, 1))

        flattened_directed_shps.append(cl)

    for ind, i in enumerate(flattened_directed_shps):
        curr_pw = get_path_weight(G, i)
        fname = "%s_scaffold_path_%d_aln.txt" % (outname, ind + 1)
        write_path_alignment(G, i, fname, curr_pw)

    write_path_cycles(G, flattened_directed_shps, outname + "_scaffold_paths.txt")

    # graph to cytoscape js file
    graph_dict = graphs_to_cytoscapejs_dict(G)
    fname = outname + "_data.json"
    with open(fname, 'w') as outfile:
        json.dump(graph_dict, outfile)

    logging.info("\nFinished")
    logging.info(str(datetime.datetime.now()) + "\n")
    logging.shutdown()
