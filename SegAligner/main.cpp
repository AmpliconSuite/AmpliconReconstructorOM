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

using namespace std;

//const float fp_p = 5000.0;
int lookback = 5;
int min_map_len = 4;
bool limit_lookback = true;
bool swap_b_x = false;

/**
 * Opens and parses a CMAP text file, stores cmap as pair which is given as input
 * @param fname
 * @param cmap_pair
 */
void parse_cmap(const string &fname, map<int,vector<float>> &cmap_pair) {
    vector<float> *new_map = new vector<float>;
    string line;
    vector<string> tokens;
    ifstream infile (fname);
    size_t pos = 0;
    string token;
    int curr_id = 0;
    float curr_pos = 0;
    if (infile.is_open()) {
        while (getline (infile,line)) {
            int counter = 0;
            if (line[0] == '#') {
                continue;
            }
            while ((pos = line.find('\t')) != string::npos) {
                token = line.substr(0,pos);
                if (counter == 0) {
                    if (curr_id != stoi(token)) {
                        curr_id = stoi(token);
                        cmap_pair[curr_id] = *new_map;
                    }
                } else if (counter == 5) {
                    curr_pos = stof(token);
                    break;
                }
                line.erase(0,pos+1);
                counter++;
            }
            cmap_pair[curr_id].push_back(curr_pos);
        }
        infile.close();
    }
}

///**
// * Prints a cmap map by map-id to cout
// * @param cmap_map
// */
//void cmap_map_to_string(map<string,vector<float>> cmap_map) {
//    for (map<string,vector<float>>::const_iterator iter = cmap_map.begin(); iter != cmap_map.end(); ++iter) {
//        cout << "CMAP ID " << iter->first << endl;
//        vector<float> curr_posns = iter->second;
//        for (vector<float>::const_iterator iter2 = curr_posns.begin(); iter2 != curr_posns.end(); ++iter2) {
//            cout << "POSN " << *iter2 << endl;
//        }
//    }
//}

/**
 * Makes reverse cmap_map entries and adds them into the map
 * @param cmap_map
 */
map<int,vector<float>> make_reverse_cmap(map<int,vector<float>> &cmap_map, int min_map_len) {
    map<int,vector<float>> full_cmap_map;
    vector<float> *new_map;
    int map_id;
    int rev_map_id;
    float map_length;
    vector<float> curr_posns;
    for (auto i = cmap_map.begin(); i != cmap_map.end(); ++i) {
        curr_posns = i->second;
        map_id = i->first;
        //full_cmap_map[map_id] = i->second;
        map_length = curr_posns.size();
        //to do: check if length of curr_posns > 1
        if (map_length <= min_map_len + 1) {
            continue;
        }
        rev_map_id = -1*map_id;
        new_map = new vector<float>;
        full_cmap_map[map_id] = curr_posns;
        full_cmap_map[rev_map_id] = *new_map;
        for (int j = curr_posns.size() - 2; j > -1; --j) {
            full_cmap_map[rev_map_id].push_back(map_length - curr_posns[j]);
        }
        full_cmap_map[rev_map_id].push_back(map_length);
    }
    return full_cmap_map;
}

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
    float delta = powf(fabs((b_posns[j_ind] - b_posns[i_ind]) - (x_posns[q_ind] - x_posns[p_ind])),1.2f);
    float fn_t = 5000.0f*(j_ind - (i_ind+1));
    float fp_t = 5000.0f*exp_x_labels;
    return 10000.0f - (fn_t + fp_t + delta);
}


void dp_align(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous, vector<float> &b_posns,
              vector<float> &x_posns, int x, int lookback, set<array<int,3>> &used_pairings,
              vector<float> &x_collapse_probs, vector<float> &b_collapse_probs) {

//    array<int, 3> prev;
//    array<int, 2> empty = {-1, -1};
    float newScore;
    int i_start;
    int p_start;

    for (int j_ind = 0; j_ind < b_posns.size() - 1; j_ind++) {
        i_start = max(0, j_ind - lookback);
        for (int q_ind = 0; q_ind < x_posns.size() - 1; ++q_ind) {
            previous[j_ind][q_ind] = {-1, -1};

            if (used_pairings.find({x,j_ind,q_ind}) != used_pairings.end()) {
                continue;
            }

            p_start = max(0, q_ind - lookback);
            for (int i_ind = i_start; i_ind < j_ind; i_ind++) {
                for (int p_ind = p_start; p_ind < q_ind; p_ind++) {
//                    if (used_pairings.find({x, i_ind, p_ind}) == used_pairings.end()) {
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
//    cout << best_pair[0] << " " << best_pair[1] << " " << best_pair[2] <<  " " << best_score << "\n";
    return best_pair;
}

void zero_out_aln(vector<vector<float>> &S,  vector<vector<array<int,2>>> &previous, int best_j, int best_x) {
    array<int,2> new_inds = previous[best_j][best_x];
    S[best_j][best_x] = -numeric_limits<float>::infinity();
    while (new_inds[0] != -1) {
        S[new_inds[0]][new_inds[1]] = -numeric_limits<float>::infinity();
        new_inds = previous[new_inds[0]][new_inds[1]];
    }
}

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
    bool hit;
    float curr_score;
    for (int i = 0; i < n; i++) {
        hit = false;
        while (!hit && (last_high_ind > -1)) {
            last_high_ind-=1;
            curr_score = S[get<1>(all_scores[last_high_ind])][get<2>(all_scores[last_high_ind])];
            if (curr_score > 0) {
                hit = true;
                zero_out_aln(S,previous,get<1>(all_scores[last_high_ind]),get<2>(all_scores[last_high_ind]));
                scoring_vector.emplace_back(x, contig_id, curr_score);
            }
        }
    }
}

/**
 * Open an output file and print alignment path there
 * @param S
 * @param previous
 * @param start
 * @param b_cid
 */
void print_alignment(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous, array<int,2> &start,
                     int contig_id, int x, map<int,vector<float>> cmaps_ref, float &evalue,
                     string &outname, set<array<int,3>> &used_pairings) {

    ofstream outfile;
    int seg_id = x;
    char direction;
    if (seg_id < 0) {
        direction = '-';
    } else {
        direction = '+';
    }
    outfile.open(outname);
    outfile << "#seg_seq\ttotal_score\tcircular\n";
    outfile << "#"  << abs(x) << direction << "\t" << S[start[0]][start[1]] << "\tFalse\n";
    outfile << "#contig_id\tseg_id\tcontig_label\tseg_label\tcontig_dir\tseg_dir\tseg_aln_number\tscore\tscore_delta\n";
    //b_label, seg label, score
    vector<tuple<int,int,float>> aln_list;

    //backtrack and build up list
    array<int,2> curr = start;
    while (curr[0] != -1) {
//        aln_list.emplace_back(make_tuple(curr[0],curr[1],S[curr[0]][curr[1]]));
        aln_list.emplace_back(curr[0],curr[1],S[curr[0]][curr[1]]);
        used_pairings.insert({x,curr[0],curr[1]});
        curr = previous[curr[0]][curr[1]];
    }

    int n_labs;
    float prev_score = 0.0;
    float curr_score;
    float score_delta;
    int curr_lab;
    for (int i = aln_list.size()-1; i > -1; --i) {
        curr_score = get<2>(aln_list[i]);
        score_delta = curr_score - prev_score;
        curr_lab = get<1>(aln_list[i]);
        if (seg_id < 0) {
//            direction = '-';
            //translate the reverse label number
//            seg_id = -1*seg_id;
            n_labs = cmaps_ref[seg_id].size();
            curr_lab = n_labs - curr_lab - 1;

        } else {
            curr_lab+=1;
        }
//        fprintf(outfile,"%d\t%d\t%d\t%d\t+\t%c\t0\t%f\t%f\n",contig_id,seg_id,curr_lab,get<0>(aln_list[i]),direction,curr_score,score_delta);
        outfile << contig_id << "\t" << abs(seg_id) << "\t" << get<0>(aln_list[i])+1 << "\t" << curr_lab << "\t+\t"
                << direction << "\t0\t" << curr_score << "\t" << score_delta << "\n";

        prev_score = curr_score;
    }

    outfile << flush;
    outfile.close();
}

void init_semiglobal(vector<vector<float>> &S, int x_size, int b_size, int x) {
    for (int q_ind = 0; q_ind < x_size - 1; q_ind++) {
        S[0][q_ind] = 0.0;
    }

    for (int j_ind = 0; j_ind < b_size - 1; j_ind++) {
        S[j_ind][0] = 0.0;
    }

    for (int j_ind = 1; j_ind < b_size - 1; j_ind++) {
        for (int q_ind = 1; q_ind < x_size - 1; ++q_ind) {
            S[j_ind][q_ind] = -numeric_limits<float>::infinity();
        }
    }
}

void init_local_aln(vector<vector<float>> &S, int x_size, int b_size, int x) {
    for (int j_ind = 0; j_ind < b_size - 1; j_ind++) {
        for (int q_ind = 0; q_ind < x_size - 1; ++q_ind) {
            S[j_ind][q_ind] = 0;
        }
    }
}


//Parse the score thresholds
map<int,float>parse_score_thresholds(string &score_thresh_file) {
    map<int,float>score_thresholds;
    ifstream infile(score_thresh_file);
    int seg_id;
    float score;
    while (infile >> seg_id >> score) {
        if (score > 0) {
//            cout << seg_id << "\n";
            score_thresholds[seg_id] = score;
            score_thresholds[-1 * seg_id] = score;
        }
    }
    return score_thresholds;
}

//Parse list of contigs to aln with
set<int>parse_contig_list(string &contig_list_file) {
    set<int>contig_set;
    ifstream infile(contig_list_file);
    int contig_id;
    while (infile >> contig_id) {
        contig_set.insert(contig_id);
    }
    return contig_set;
}

/*
 *MAIN
 */
int main (int argc, char *argv[]) {
    if (argc < 3) {
        cout << "wrong number of arguments" << endl;
        exit(1);
    }

    //timing
    clock_t begin = clock();
    double elapsed_secs;
    clock_t end;

    //positional args
    string ref_cmap_file = argv[1];
    string contig_cmap_file = argv[2];
    min_map_len = atoi(argv[3]);

    //non-positional arg inits:
    bool scoring_mode = false;
    bool local_aln = false;
    bool detection = false;
    string sample_prefix;
    string score_thresh_file;
    string contig_list_file;

    map<int,float>score_thresholds;
    vector<tuple<int,int,float>>scoring_vector;
    set<array<int,3>>used_pairings;
    int total_labels = 0;
    string detection_label = "";
    string flipped = "";

    //make segs cmap
    map<int,vector<float>> cmaps_ref_raw;
    ;
    parse_cmap(ref_cmap_file,cmaps_ref_raw);
    //    cmap_map_to_string(cmap_map) //DEBUG
    //add reverse segs
    map<int,vector<float>> cmaps_ref = make_reverse_cmap(cmaps_ref_raw, min_map_len);

    //make contigs cmaps
    map<int,vector<float>> cmaps_contig;
    parse_cmap(contig_cmap_file, cmaps_contig);

    //find collapse probs
    map<int,vector<float>> non_collapse_prob_map = non_collapse_probs(cmaps_ref);
    map<int,vector<float>> non_collapse_prob_map_contig = non_collapse_probs(cmaps_contig);

    for (int i = 3; i < argc; ++i) {
        if (string(argv[i]) == "-nl") {
            limit_lookback = false;
        } else if (string(argv[i]) == "-scoring") {
            scoring_mode = true;
            cout << "Scoring mode ON" << endl;
        } else if (string(argv[i]) == "-detection") {
                local_aln = true;
                detection = true;
                swap_b_x = true;
                detection_label = "detection_";
                cout << "detection mode ON (local ON)" << endl;
        } else if (string(argv[i]) == "-alnref") {
            local_aln = true;
            swap_b_x = true;
            flipped = "_flipped";
            cout << "detection mode ON (local ON)" << endl;
        } else if (string(argv[i]).rfind("-score_thresh=",0) == 0) {
            score_thresh_file = string(argv[i]).substr(string(argv[i]).find('=')+1);
        } else if(string(argv[i]).rfind("-prefix=",0) == 0) {
            sample_prefix = string(argv[i]).substr(string(argv[i]).find('=')+1);
        } else if(string(argv[i]).rfind("-contig_list=",0) == 0) {
            contig_list_file = string(argv[i]).substr(string(argv[i]).find('=')+1);
        } else if (string(argv[i]) == "-local") {
            local_aln = true;
            cout << "local alignment mode" << endl;
        }
    }
    if (!limit_lookback) {
        lookback = 9999999;
        cout << "lookback limit OFF" << endl;
    }

    //Get Score thresholds
    if (!score_thresh_file.empty()) {
        score_thresholds = parse_score_thresholds(score_thresh_file);
    } else {
        for (auto x = cmaps_ref.begin(); x != cmaps_ref.end(); ++x) {
            for (int ncp = 0; ncp < non_collapse_prob_map[x->first].size(); ncp++) {
                score_thresholds[x->first]+=(8000.0f * non_collapse_prob_map[x->first][ncp]);
            }
        //cout << x->first << " " << score_thresholds[x->first] << "\n";
        }
    }

    //Get contigs to align with
    set<int> contig_set;
    if (!contig_list_file.empty()) {
        contig_set = parse_contig_list(contig_list_file);
    } else {
        for (auto key = cmaps_contig.begin(); key != cmaps_contig.end(); ++key) {
            contig_set.insert(key->first);
        }
    }

    //create alignment data structures
//    map<array<int,3>,float> S;
    vector<vector<float>> S;
//    map<array<int,2>,array<int,2>> previous;
    vector<vector<array<int,2>>> previous;
    vector<float> curr_ncp_vector;
    vector<float> contig_ncp_vector;
    array<int, 2> best_start = {-1,-1};
    for (auto b = cmaps_contig.begin(); b != cmaps_contig.end(); ++b) {
        total_labels+=b->second.size();
    }
    cout << "total labels in contig set: " << total_labels << endl;

    //create variables for the segment and contig iterations
    int contig_id;
    vector<float> contig_posns;
    int x;
    vector<float> x_posns;
    float curr_best_score;
    int contig_aligned_count;
    string outname = sample_prefix;
    string seg_id;
    int alns_to_get;
    cout << "Running SegAligner on " << cmaps_ref.size() << " segments and " << cmaps_contig.size() << " contigs\n";
    for (auto x_map = cmaps_ref.begin(); x_map != cmaps_ref.end(); ++x_map) {
        x = x_map->first;
        x_posns = x_map->second;
//        cout << x << endl;
        contig_aligned_count = 0;
        curr_ncp_vector = non_collapse_prob_map[x];
        for (auto contig_map = cmaps_contig.begin(); contig_map != cmaps_contig.end();) {
            contig_id = contig_map->first;
            contig_posns = contig_map->second;
            if (contig_set.find(contig_id) == contig_set.end()) {
                ++contig_map;
                continue;
            }

            //initialize alignment data structures
            for (int i = 0; i < contig_posns.size()-1; i++) {
//                S.push_back(vector<float>(x_posns.size()-1));
//                previous.push_back(vector<array<int,2>>(x_posns.size()-1));
                S.emplace_back(x_posns.size()-1);
                previous.emplace_back(x_posns.size()-1);
            }

            if (local_aln) {
                init_local_aln(S,x_posns.size(),contig_posns.size(),x);
            } else {
                init_semiglobal(S,x_posns.size(),contig_posns.size(),x);
            }
            contig_ncp_vector = non_collapse_prob_map_contig[contig_id];

            dp_align(S, previous, contig_posns, x_posns, x, lookback, used_pairings, curr_ncp_vector, contig_ncp_vector);
            best_start = get_backtrack_start(S, contig_posns.size(), x_posns.size());
            curr_best_score = S[best_start[0]][best_start[1]];

            //if mode "scoring", write this score into a list somewhere
            if (scoring_mode) {
                scoring_vector.emplace_back(x, contig_id, curr_best_score);
                ++contig_map;

            } else if (detection) {
                //figure out how many alns to get for this thing
                alns_to_get =  int(roundf(1000.0f * float(contig_posns.size())/float(total_labels))) + 1;
//                cout << contig_id << " " << x << " contig has len: " << contig_posns.size() << " getting " << alns_to_get << endl;
                multi_backtrack_scores(S,previous,x,contig_id,contig_posns.size(),x_posns.size(),alns_to_get,scoring_vector);
                ++contig_map;

            } else {
                //if it passes the e-value write it, and do not update the contig
                if (curr_best_score >= score_thresholds[x]) {
                    contig_aligned_count+=1;
                    seg_id = to_string(abs(x));
                    if (x < 0){
                        seg_id+="_r";
                    }
                    outname = sample_prefix + "segalign_" + to_string(contig_id) + "_" + seg_id + "_" + to_string(contig_aligned_count) + flipped + "_aln.txt";
                    print_alignment(S,previous,best_start,contig_id,x,cmaps_ref,score_thresholds[x],outname,used_pairings);
//                    --contig_map;
                    if ((contig_aligned_count > 4) && local_aln) {
                        ++contig_map;
                    }

                } else {
                    used_pairings.clear();
                    ++contig_map;
                }
            }
            S.clear();
            previous.clear();
        }
    }
    if (scoring_mode || detection) {
        ofstream outfile;
        outfile.open(outname + detection_label + "scores_" + to_string(contig_id) + "_contig_scores.txt");
        outfile << "#seg\tcontig\tscore\n";
        for (int i = 0; i < scoring_vector.size(); i++) {
            outfile << get<0>(scoring_vector[i]) << "\t" << get<1>(scoring_vector[i]) << "\t"
                    << get<2>(scoring_vector[i]) << "\n";
        }
        outfile << flush;
        outfile.close();
    }

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("elapsed seconds for run %f\n",elapsed_secs);
}