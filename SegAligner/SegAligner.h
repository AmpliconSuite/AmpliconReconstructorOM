using namespace std;

map<int,float> compute_score_thresholds(vector<tuple<int,int,float>> &full_scoring_vector, map<int,vector<float>> &cmaps_segs,
                                        map<int,vector<float>> &cmaps_contigs, double pvalue);

pair<float,float> compute_partial_score_threshold(float full_threshold, const vector<tuple<int,int,float>> &aln_list,
                                                  vector<float> &contig_cmap);

map<int,vector<float>> non_collapse_probs(map<int,vector<float>> &segs_cmaps, int inst_gen);

float score_f(vector<float> &b_posns, vector<float> &x_posns, int i_ind, int j_ind, int p_ind, int q_ind,
              vector<float> &x_collapse_probs);

void dp_align(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous, vector<float> &b_posns,
              vector<float> &x_posns, bool fitting_aln, int lookback, set<int> &used_labels,
              vector<float> &x_collapse_probs, vector<float> &b_collapse_probs, bool swap_b_x);

array<int,2> local_backtrack_start(vector<vector<float>> &S,int b_len,int x_len);

array<int,2> score_backtrack_start(vector<vector<float>> &S,int b_len,int x_len);

array<int,2> sg_aln_backtrack_start(vector<vector<float>> &S,int b_len,int x_len);

void zero_out_aln(vector<vector<float>> &S,  vector<vector<array<int,2>>> &previous, int best_j, int best_x);

void multi_backtrack_scores(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
                            int x, int contig_id, int b_len, int x_len, int n, vector<tuple<int,int,float>> &scoring_vector);

bool check_reused_labels(vector<tuple<int,int,float>> &aln_list, set<int> &used_label_set);

void add_used_labels(vector<tuple<int,int,float>> &aln_list, int c_id, map<int,set<int>> &contig_used_label_map);

void init_semiglobal_aln(vector<vector<float>> &S, int x_size, int b_size, int x);

void init_local_aln(vector<vector<float>> &S, int x_size, int b_size, int x);

void init_fitting_aln(vector<vector<float>> &S, int x_size, int b_size, int x);

//void init_semiglobal_scoring(vector<vector<float>> &S, int x_size, int b_size, int x);

pair<double,double> linreg(vector<float> &x, vector<int> &y);

float compute_aln_median(vector<tuple<int, int, float>> aln_list);

float compute_aln_min(vector<tuple<int, int, float>> aln_list);