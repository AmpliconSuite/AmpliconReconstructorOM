#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
//#include <string.h>
#include <vector>
#include <map>
#include <ctime>
#include <cmath>
#include <set>
#include <limits>
#include <numeric>

#include "SAHelper.h"

using namespace std;

//const float fp_p = 5000.0;

//iterate over in-silico reference and estimate a probability that a given label will appear as collapsed with a left
// or right neighbor
map<int,vector<float>> non_collapse_probs(map<int,vector<float>> &segs_cmaps) {
    map<int,vector<tuple<float,float>>> non_collapse_pairs;
    map<int,vector<float>> non_collapse_prob_map;
    int id;
    int y_len;
    float dist;
    float non_collapse = 1.0;
    vector<float> y_posns;
    for (auto y = segs_cmaps.begin(); y != segs_cmaps.end(); ++y) {
        id = y->first;
        y_posns = y->second;
        y_len = int(y->second.size()) - 1;
        non_collapse_pairs[id].emplace_back(1.0,1.0);
        for (int i = 1; i < y_len; ++i) {
            non_collapse_pairs[id].emplace_back(1.0,1.0);
            dist = y_posns[i] - y_posns[i-1];
            non_collapse = fmin(1.0f,powf(dist,4)/powf(2000,4));
            non_collapse_pairs[id][i] = make_pair(non_collapse,1.0);
            non_collapse_pairs[id][i-1] = make_pair(get<0>(non_collapse_pairs[id][i-1]),non_collapse);
            non_collapse_prob_map[id].push_back(get<0>(non_collapse_pairs[id][i-1])*non_collapse);
        }
        non_collapse_prob_map[id].push_back(get<0>(non_collapse_pairs[id][y_len-1])*non_collapse);
    }
    return non_collapse_prob_map;
}


/**
 * scores a matching region
 * @param b_posns
 * @param x_posns
 * @param x
 * @param i_ind
 * @param j_ind
 * @param p_ind
 * @param q_ind
 * @param non_collapse_prob_map
 * @return
 */
float score_f(vector<float> &b_posns, vector<float> &x_posns, int i_ind, int j_ind, int p_ind, int q_ind,
              vector<float> &x_collapse_probs) {

    float exp_x_labels = accumulate(x_collapse_probs.begin() + p_ind + 1, x_collapse_probs.begin() + q_ind, 0.0f);
    float discrep = fabs((b_posns[j_ind] - b_posns[i_ind]) - (x_posns[q_ind] - x_posns[p_ind]));
    float delta = powf(discrep,1.2f);
    float fn_t = 5000.0f*(j_ind - (i_ind+1));
    float fp_t = 5000.0f*exp_x_labels;
    return 10000.0f - (fn_t + fp_t + delta);
}


void dp_align(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous, vector<float> &b_posns,
              vector<float> &x_posns, int x, int lookback, set<int> &used_labels,
              vector<float> &x_collapse_probs, vector<float> &b_collapse_probs, bool swap_b_x) {

    for (int j_ind = 0; j_ind < b_posns.size() - 1; j_ind++) {
        int i_start = max(0, j_ind - lookback);
        for (int q_ind = 0; q_ind < x_posns.size() - 1; ++q_ind) {
            previous[j_ind][q_ind] = {-1, -1};

            if (used_labels.find(j_ind) != used_labels.end()) {
                continue;
            }

            int p_start = max(0, q_ind - lookback);
            for (int i_ind = i_start; i_ind < j_ind; i_ind++) {
                for (int p_ind = p_start; p_ind < q_ind; p_ind++) {
                    float newScore;
                    if (!swap_b_x) {
                        newScore = S[i_ind][p_ind] +
                                   score_f(b_posns, x_posns, i_ind, j_ind, p_ind, q_ind, x_collapse_probs);
                    } else {
                        newScore = S[i_ind][p_ind] +
                                   score_f(x_posns, b_posns, p_ind, q_ind, i_ind, j_ind, b_collapse_probs);
                    }

                    if (newScore > S[j_ind][q_ind]) {
                        S[j_ind][q_ind] = newScore;
                        previous[j_ind][q_ind] = {i_ind, p_ind};
                    }

                }
            }
        }
    }
}

/**
 * Find location to begin backtracking
 * @param cmap_map_ref
 * @param S
 * @param b_last
 * @return map key for backtracking map which has best score
 */
array<int,2> get_backtrack_start(vector<vector<float>> &S,int b_len,int x_len) {
    float best_score = -numeric_limits<float>::infinity();
    array<int,2> best_pair = {-1,-1};
    for (int j_ind = 0; j_ind < b_len-1; j_ind++) {
        for (int x_ind = 0; x_ind < x_len-1; x_ind++) {
            if (S[j_ind][x_ind] > best_score) {
                best_score = S[j_ind][x_ind];
                best_pair = {j_ind,x_ind};
            }
        }
    }
    return best_pair;
}

array<int,2> score_backtrack_start(vector<vector<float>> &S,int b_len,int x_len) {
    float best_score = -numeric_limits<float>::infinity();
    array<int,2> best_pair = {0,0};

    for (int j_ind = 0; j_ind < b_len-1; j_ind++) {
        for (int x_end_ind = x_len-3; x_end_ind < x_len-2; x_end_ind++)
        if (S[j_ind][x_end_ind] > best_score) {
            best_score = S[j_ind][x_end_ind];
            best_pair = {j_ind, x_end_ind};
        }
    }
    return best_pair;
}

//get backtrack start for tip alignment case
array<int,2> tip_aln_backtrack_start(vector<vector<float>> &S,int b_len,int x_len) {
    float best_score = -numeric_limits<float>::infinity();
    array<int,2> best_pair = {-1,-1};

    //check all pairings with ends of contig
    int j_end_ind = b_len - 2;
    for (int x_ind = 0; x_ind < x_len-1; x_ind++) {
        if (S[j_end_ind][x_ind] > best_score) {
            best_score = S[j_end_ind][x_ind];
            best_pair = {j_end_ind, x_ind};
        }
    }
    int x_end_ind = x_len - 2;
    for (int j_ind = 0; j_ind < b_len-1; j_ind++) {
        if (S[j_ind][x_end_ind] > best_score) {
            best_score = S[j_ind][x_end_ind];
            best_pair = {j_ind, x_end_ind};
        }
    }

    return best_pair;
}

//remove a used alignment from the scoring matrix
void zero_out_aln(vector<vector<float>> &S,  vector<vector<array<int,2>>> &previous, int best_j, int best_x) {
    array<int,2> new_inds = previous[best_j][best_x];
    S[best_j][best_x] = -numeric_limits<float>::infinity();
    while (new_inds[0] != -1) {
        S[new_inds[0]][new_inds[1]] = -numeric_limits<float>::infinity();
        new_inds = previous[new_inds[0]][new_inds[1]];
    }
}

//get multiple backtracking scores from a single scoring matrix
void multi_backtrack_scores(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
                            int x, int contig_id, int b_len, int x_len, int n, vector<tuple<int,int,float>> &scoring_vector) {

    vector<tuple<float,int,int>> all_scores;
    for (int j_ind = 0; j_ind < b_len-1; j_ind++) {
        for (int x_ind = 0; x_ind < x_len-1; x_ind++) {
            if (S[j_ind][x_ind] > 0) {
                all_scores.emplace_back(S[j_ind][x_ind], j_ind, x_ind);
            }
        }
    }
    sort(all_scores.begin(),all_scores.end());
    int last_high_ind = all_scores.size()-1;
    for (int i = 0; i < n; i++) {
        bool hit = false;
        while (!hit && (last_high_ind > -1)) {
            last_high_ind-=1;
            float curr_score = S[get<1>(all_scores[last_high_ind])][get<2>(all_scores[last_high_ind])];
            if (curr_score > 0) {
                hit = true;
                zero_out_aln(S,previous,get<1>(all_scores[last_high_ind]),get<2>(all_scores[last_high_ind]));
                scoring_vector.emplace_back(x, contig_id, curr_score);
            }
        }
    }
}

//during tip alignment check if an alignment overlaps previously aligned part of contig
bool check_reused_labels(vector<tuple<int,int,float>> &aln_list, set<int> &used_label_set) {
    for (auto item = aln_list.begin(); item != aln_list.end(); ++item) {
        if (used_label_set.find(get<1>(*item)) != used_label_set.end()) {
            cout << "hit an overlap " << get<1>(*item) << "\n";
            return true;
        }
    }
    return false;
}

//during tip alignment add previously kept contig labels so they aren't used again
void add_used_labels(vector<tuple<int,int,float>> &aln_list, int c_id, map<int,set<int>> &contig_used_label_map) {
    for (auto item = aln_list.begin(); item != aln_list.end(); ++item) {
        contig_used_label_map[c_id].insert(get<1>(*item));
    }
}

//initialize the scoring matrix for a semiglobal alignment
void init_semiglobal(vector<vector<float>> &S, int x_size, int b_size, int x) {
    for (int q_ind = 0; q_ind < x_size - 1; q_ind++) {
        S[0][q_ind] = 0.0;
    }

    for (int j_ind = 0; j_ind < b_size - 1; j_ind++) {
        S[j_ind][0] = 0.0;
        S[j_ind][1] = 0.0;
    }

    for (int j_ind = 1; j_ind < b_size - 1; j_ind++) {
        for (int q_ind = 2; q_ind < x_size - 1; ++q_ind) {
            S[j_ind][q_ind] = -numeric_limits<float>::infinity();
        }
    }
}

//initialize the scoring matrix for a local alignment
void init_local_aln(vector<vector<float>> &S, int x_size, int b_size, int x) {
    for (int j_ind = 0; j_ind < b_size - 1; j_ind++) {
        for (int q_ind = 0; q_ind < x_size - 1; ++q_ind) {
            S[j_ind][q_ind] = 0;
        }
    }
}

/*
 * Get the slope and intercept for linear regression
 */
pair<double,double> linreg(vector<float> &x, vector<int> &y) {
    double n = x.size();
    double xbar = accumulate(x.begin(), x.end(), 0.0) / n;
    double ybar = accumulate(y.begin(), y.end(), 0.0) / n;
    double cov_xy = 0.0;
    double var_x = 0.0;
    for (int i=0; i<n; i++) {
        cov_xy+= (x[i] - xbar) * (y[i] - ybar);
        var_x+= (x[i] - xbar) * (x[i] - xbar);
    }

    double beta = cov_xy/var_x;
    double intercept = ybar - beta*xbar;
    return make_pair(beta,intercept);

}

float compute_aln_median(vector<tuple<int, int, float>> aln_list) {
    //get score median
    vector<float> score_deltas;
    float prev = 0.0;
    for (int i = aln_list.size()-2; i > -1; i--) {
        float diff = get<2>(aln_list[i]) - prev;
        score_deltas.push_back(diff);
        prev = get<2>(aln_list[i]);
    }
    sort(score_deltas.begin(),score_deltas.end());
    float score_median;
    int mid;
    if (score_deltas.size() % 2) {
        mid = int(score_deltas.size())/2;
        score_median = score_deltas[mid];

    } else {
        mid = int(score_deltas.size())/2;
        score_median = (score_deltas[mid] + score_deltas[mid-1])/2;
    }

    return score_median;
}