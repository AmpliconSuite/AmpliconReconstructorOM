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

void parse_cmap(const string &fname, map<int,vector<float>> &cmap_pair) {
    auto *new_map = new vector<float>;
    string line;
    ifstream infile(fname);
    size_t pos = 0;
    int curr_id = 0;
    float curr_pos = 0;
    if (infile.is_open()) {
        while (getline (infile,line)) {
            int counter = 0;
            if (line[0] == '#') {
                continue;
            }
            while ((pos = line.find('\t')) != string::npos) {
                string token = line.substr(0,pos);
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

void cmap_map_to_string(map<string,vector<float>> cmap_map) {
    for (map<string,vector<float>>::const_iterator iter = cmap_map.begin(); iter != cmap_map.end(); ++iter) {
        cout << "CMAP ID " << iter->first << endl;
        vector<float> curr_posns = iter->second;
        for (vector<float>::const_iterator iter2 = curr_posns.begin(); iter2 != curr_posns.end(); ++iter2) {
            cout << "POSN " << *iter2 << endl;
        }
    }
}

map<int,vector<float>> make_reverse_cmap(map<int,vector<float>> &cmap_map, int min_map_len) {
    map<int,vector<float>> full_cmap_map;
    for (auto i: cmap_map) {
        vector<float> curr_posns = i.second;
        int map_id = i.first;
        //full_cmap_map[map_id] = i->second;
        float map_length = curr_posns.size();
        //to do: check if length of curr_posns > 1
        if (map_length <= min_map_len + 1) {
            continue;
        }
        int rev_map_id = -1*map_id;
        full_cmap_map[map_id] = curr_posns;
        full_cmap_map[rev_map_id] = vector<float>();
        for (int j = curr_posns.size() - 2; j > -1; --j) {
            full_cmap_map[rev_map_id].push_back(map_length - curr_posns[j]);
        }
        full_cmap_map[rev_map_id].push_back(map_length);
    }
    return full_cmap_map;
}

vector<tuple<int,int,float>> get_aln_list(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
        array<int,2> &start) {

    //contig_label, segment_label, score
    vector<tuple<int,int,float>> aln_list;
    //backtrack and build up list
    array<int,2> curr = start;
    while (curr[0] != -1) {
        aln_list.emplace_back(curr[0],curr[1],S[curr[0]][curr[1]]);
        curr = previous[curr[0]][curr[1]];
    }
    return aln_list;
}

void print_alignment(vector<vector<float>> &S, vector<vector<array<int,2>>> &previous,
        vector<tuple<int,int,float>> &aln_list, int contig_id, int seg_id, map<int,vector<float>> cmaps_ref,
        string &outname, set<array<int,3>> &used_pairings, float curr_best_score) {

    ofstream outfile;
    char direction;
    if (seg_id < 0) {
        direction = '-';
    } else {
        direction = '+';
    }
    outfile.open(outname);
    outfile << "#seg_seq\ttotal_score\tcircular\n";
    outfile << "#"  << abs(seg_id) << direction << "\t" << curr_best_score << "\tFalse\n";
    outfile << "#contig_id\tseg_id\tcontig_label\tseg_label\tcontig_dir\tseg_dir\tseg_aln_number\tscore\tscore_delta\n";

    float prev_score = 0.0;
    for (int i = aln_list.size()-1; i > -1; --i) {
        used_pairings.insert({seg_id,get<0>(aln_list[i]),get<1>(aln_list[i])});
        float curr_score = get<2>(aln_list[i]);
        float score_delta = curr_score - prev_score;
        int curr_lab = get<1>(aln_list[i]);
        if (seg_id < 0) {
            //translate the reverse label number
            int n_labs = cmaps_ref[seg_id].size();
            curr_lab = n_labs - curr_lab - 1;

        } else {
            curr_lab+=1;
        }
        outfile << contig_id << "\t" << abs(seg_id) << "\t" << get<0>(aln_list[i])+1 << "\t" << curr_lab << "\t+\t"
                << direction << "\t0\t" << curr_score << "\t" << score_delta << "\n";

        prev_score = curr_score;
    }

    outfile << flush;
    outfile.close();
}

//map<int,float>parse_score_thresholds(string &score_thresh_file) {
//    map<int,float>score_thresholds;
//    ifstream infile(score_thresh_file);
//    int seg_id;
//    float score;
//    while (infile >> seg_id >> score) {
//        if (score > 0) {
//            score_thresholds[seg_id] = score;
//            score_thresholds[-1 * seg_id] = score;
//        }
//    }
//    return score_thresholds;
//}

set<int>parse_contig_list(string &contig_list_file) {
    set<int>contig_set;
    ifstream infile(contig_list_file);
    int contig_id;
    while (infile >> contig_id) {
        contig_set.insert(contig_id);
    }
    return contig_set;
}

void parse_used_labels(string &used_labels_file, map<int,set<int>> &contig_used_label_map) {
    ifstream infile(used_labels_file);
    string line;
    size_t pos;
    int contig_id;
    while (getline(infile,line)) {
        if (line[0] == '#') {
            continue;
        }
        int counter = 0;
        while ((pos = line.find(' ')) != string::npos) {
            string token = line.substr(0, pos);
            if (counter == 0) {
                contig_id = stoi(token);
            } else {
                contig_used_label_map[contig_id].insert(stoi(token));
            }

            line.erase(0, pos + 1);
            counter++;
        }

    }
    infile.close();
}

void write_score_thresholds(map<int,float> score_thresholds, string full_prefix) {
    //write scoring distribution
    cout << "Writing scoring distribution.\n";
    ofstream outfile;
    outfile.open(full_prefix + "_score_thresholds.txt");
    outfile << "#seg\tscore\n";
    for (auto &e: score_thresholds) {
        outfile << e.first << "\t" << e.second << "\n";
    }
    outfile << flush;
    outfile.close();
}

bool null_intersection(vector<tuple<int,int,float>> &aln_list, set<int> &s) {
    for (auto &e: aln_list) {
        if (s.find(get<0>(e)) != s.end()){
            return false;
        }
    }
    return true;
}

void write_all_scores(vector<tuple<int,int,float>> &full_scores, string &prefix) {
    ofstream outfile;
    outfile.open(prefix + "_all_scores.txt");
    outfile << "#seg\tcontig\tscore\n";
    for (auto &e: full_scores) {
        outfile << get<0>(e) << "\t" << get<1>(e) << "\t" << get<2>(e) << "\n";
    }
    outfile << flush;
    outfile.close();
}