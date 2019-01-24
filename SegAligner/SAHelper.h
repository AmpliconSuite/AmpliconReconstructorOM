using namespace std;

/**
 * Opens and parses a CMAP text file, stores cmap as pair which is given as input
 * @param fname
 * @param cmap_pair
 */
void parse_cmap(const string &fname, map<int,vector<float>> &cmap_pair);

///**
// * Prints a cmap map by map-id to cout
// */
void cmap_map_to_string(map<int,vector<float>> cmap_map);


//Makes reverse cmap_map entries and adds them into the map
map<int,vector<float>> make_reverse_cmap(map<int,vector<float>> &cmap_map, int min_map_len);

//get the alignment and mark labels as used
vector<tuple<int,int,float>> get_aln_list(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
        array<int,2> &start);

//Open an output file and print alignment path there, mark the output alignments as used
void print_alignment(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
        vector<tuple<int,int,float>> &aln_list, int contig_id, int seg_id, map<int,vector<float>> cmaps_ref,
        string &outname, map<int,set<int>> &contig_used_label_map, float curr_best_score);

////Parse the score thresholds
//map<int,float> parse_score_thresholds(string &score_thresh_file);

//Parse list of contigs to aln with
set<int> parse_contig_list(string &contig_list_file);

//Parse set of labels per contig which have already been paired
void parse_used_labels(string &used_labels_file, map<int,set<int>> &contig_used_label_map);

void write_score_thresholds(map<int,float> score_thresholds, string full_prefix);

//bool null_intersection(vector<tuple<int,int,float>> &aln_list, set<int> &s);

void write_all_scores(vector<tuple<int,int,float>> &full_scores, string &full_prefix);