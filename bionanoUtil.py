def median(L):
    if L:
        L = sorted(L)
        n = len(L)
        m = n - 1
        return (L[n/2] + L[m/2]) / 2.0

    return None

#parse cmap
def parse_cmap(cmapf,keep_length = False):
    cmaps = {}
    #contigCovs = {}

    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                if fD["CMapId"] not in cmaps:
                    cmaps[fD["CMapId"]] = {}
                    #contigCovs[fD["CMapId"]] = {}

                #this is not a good way to parse label channel means color channel
                if fD["LabelChannel"] == "1":
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])
                    #contigCovs[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Coverage"])
                    
                elif fD["LabelChannel"] == "0" and keep_length:
                    cmaps[fD["CMapId"]][int(fD["SiteID"])] = float(fD["Position"])

    return cmaps

# #parse the BNX file and return a dictionary of molecule ID to the molecule entry
# def parseBNX(bnxF):
#     moleculeD = {}
#     with open(bnxF) as infile:
#         for line in infile:
#             if not line.startswith("#"):
#                 fields = line.rstrip().rsplit("\t")
#                 if line.startswith("0"):
#                     currKey = fields[1]
#                 else:
#                     moleculeD[currKey] = fields[1:]

#     return moleculeD

def get_cmap_lens(cmapf):
    cmap_lens = {}
    #contigCovs = {}
    with open(cmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                if fD["CMapId"] not in cmap_lens:
                    cmap_lens[fD["CMapId"]] = float(fD["ContigLength"])

    #print cmaps
    return cmap_lens


#parse xmap
def parse_xmap(xmapf):
    detailFields = ["QryContigID","RefContigID","Orientation","Confidence","QryLen","RefLen",
    "QryStartPos","QryEndPos","RefStartPos","RefEndPos","Alignment"]
    #xmapAln = {}
    xmapPair = {}
    with open(xmapf) as infile:
        for line in infile:
            if line.startswith("#h"):
                head = line.rstrip().rsplit()[1:]

            elif not line.startswith("#"):
                fields = line.rstrip().rsplit()
                fD = dict(zip(head,fields))
                alnstring = ")" + fD["Alignment"] + "("
                #xmapAln[fD["XmapEntryID"]] = alnstring.rsplit(")(")[1:-1]

                xmapPair[fD["XmapEntryID"]] = {x:fD[x] for x in detailFields}

    return xmapPair

#read in key file 
def parse_keyfile(keyF_name):
    keyCompD = {}
    with open(keyF_name) as infile:
        for line in infile:
            if not line.startswith("#"):
                if line.startswith("CompntId"):
                    head = line.rstrip().split()
                else:
                    fields = line.rstrip().split()
                    keyCompD[fields[0]] = (fields[1],float(fields[2]))
                    
    return keyCompD

#make cmap dictionary into 0-indexed vector
def vectorize_cmaps(cmap_d):
    vectorized_dict = {}
    for y in cmap_d:
        y_posns = [cmap_d[y][k] for k in sorted(cmap_d[y].keys())]
        vectorized_dict[y] = y_posns
            
    return vectorized_dict

#take label number from a reversed segment and write it back to forward
#ASSUMES 0-BASED INDEX
def translate_reversed_label(cmaps,rev_seg_id,r_lab):
    seg_id = rev_seg_id.rsplit("_")[0]
    nlabs = len(cmaps[seg_id])
    return nlabs - r_lab

#cmaps input must include the segment artificial end label to work properly
def add_full_reverse_cmaps(cmaps,key_dict):
    #make a reverse keyfile
    iter_keys = cmaps.keys()
    for i in iter_keys:
        tot_labs = len(cmaps[i])-1
        cmap_len = cmaps[i][tot_labs+1]
        new_ID = i + "_r"
        #add new entry to key dict
        seg_rep = "|".join(key_dict[i][0].rsplit("|")[::-1])
        key_dict[new_ID] = (seg_rep,cmap_len)
        #add new entry to cmaps
        cmaps[new_ID] = {}
        for j in range(1,tot_labs+1):
            cmaps[new_ID][tot_labs - j + 1] = cmap_len - cmaps[i][j]
        
        cmaps[new_ID][tot_labs+1] = cmap_len

#parses the output from SegAligner
def parse_seg_alignment_file(alignfile):
    alignment = []
    tip_aln = True if "tip_seg_aln" in alignfile else False
    with open(alignfile) as infile:
        meta_head = infile.next().rstrip()[1:].rsplit()
        meta_vals = infile.next().rstrip()[1:].rsplit()
        meta_dict = dict(zip(meta_head,meta_vals))
        aln_head = infile.next().rstrip()[1:].rsplit()
        for line in infile:
            fields = line.rstrip().rsplit()
            alignment.append(dict(zip(aln_head,fields)))

        seg_id = meta_dict["seg_seq"][:-1]
        strand = meta_dict["seg_seq"][-1]
        tot_score = float(meta_dict["total_score"])
        seg_start = int(alignment[0]["seg_label"])
        seg_end = int(alignment[-1]["seg_label"])
        seg_ends = (seg_start,seg_end)
        contig_start = int(alignment[0]["contig_label"])
        contig_end = int(alignment[-1]["contig_label"])
        contig_ends = (contig_start,contig_end)

    return alignment[0]["contig_id"],[seg_id,seg_ends,contig_ends,strand,tot_score,alignment,tip_aln]
