#include <iostream>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <ctime>
#include <cmath>
#include <set>
#include <limits>
#include <numeric>
#include <future>

#include "SAHelper.h"
#include "SegAligner.h"

//non-positional arg instantiation:
//bool scoring_mode = false;
bool local_aln = false;
bool detection = false;
bool fitting_aln = false;
bool tip_aln = false;
bool swap_b_x = false;
int n_detect_scores = 500;
int n_threads = 1;
int min_map_len = 10;
int lookback = 6;
int inst_gen = 2;
int max_seg_contig_alns = 6;
float p_val = 0.0001;
float p_val_tip = 0.00001;
float p_val_RG =  0.000000005;
string sample_prefix = "SA_output";

//bookkeeping variable instantiation
string detection_label;
string flipped;

/*
 * method to compute the e-values
 * assumes scores are sorted
 */
map<int,float> compute_score_thresholds(vector<tuple<int,int,float>> &full_scoring_vector, map<int,vector<float>> &cmaps_segs,
        map<int,vector<float>> &cmaps_contigs, double pvalue) {

    //Organize the segment scores
    map<int,vector<float>> seg_scores;
    for (auto e: cmaps_segs) {
        seg_scores[e.first] = vector<float>();
    }

    for (auto &e: full_scoring_vector) {
        seg_scores[get<0>(e)].push_back(get<2>(e));
    }

    //sort the scores
    for (auto &e: seg_scores) {
        vector<float> curr_scores = e.second;
        sort(e.second.begin(),e.second.end());
    }

//    cout << "Finished sorting scores.\n";
//    cout << "Computing scoring threshold\n";
    map<int,float> score_thresholds;
    int right_cutoff = 25;
    int left_cuttoff_from_right = int(round(0.15*fmin(cmaps_contigs.size(),n_detect_scores)));
    double E_cutoff = -log(1 - pvalue);

    int m = 0;
    for (auto i: cmaps_contigs) {
        m+=(i.second.size()-1);
    }

    for (auto i: seg_scores) {
        int seg_id = i.first;
        int n = cmaps_segs[seg_id].size()-1;
        vector<float>::const_iterator first;
        first = i.second.begin() + (i.second.size() - (right_cutoff + left_cuttoff_from_right));
        vector<float>::const_iterator last = i.second.begin() + (i.second.size() - 1 - right_cutoff);
        vector<float> relevant_scores(first,last);
        vector<int> e_vals(relevant_scores.size());
        iota(e_vals.begin(),e_vals.end(), 1);
        reverse(e_vals.begin(),e_vals.end());

        //log the e values
        transform(e_vals.begin(),e_vals.end(),e_vals.begin(),log<int>);

        //linear regression
        double slope, intercept;
        tie(slope,intercept) = linreg(relevant_scores,e_vals);

        //compute K
        double K = exp(intercept - log(m*n));

        //compute cutoff
        double S_cutoff = -log(E_cutoff/(K*m*n))/(-slope);
        score_thresholds[seg_id] = float(S_cutoff);

    }
    return score_thresholds;
}

/*
 * Method to get the scores
 * takes a subvector
 * returns a vector
 */
vector<tuple<int,int,float>> run_SA_score(map<int,vector<float>> cmaps_segs, map<int,vector<float>> cmaps_contigs,
        map<int,vector<float>> non_collapse_prob_map, map<int,vector<float>> non_collapse_prob_map_contig,
        set<int> contig_set) {

    //variables for scoring, etc.
    vector<tuple<int,int,float>>scoring_vector;
    int total_labels = 0;
    for (auto b: cmaps_contigs) {
        total_labels+=b.second.size();
    }

    //create variables for the segment and contig iterations
    for (auto x_map: cmaps_segs) {
        int x = x_map.first;
        vector<float> x_posns = x_map.second;
        vector<float> seg_ncp_vector = non_collapse_prob_map[x];
        for (auto &contig_map: cmaps_contigs) {
            int contig_id = contig_map.first;
            vector<float> contig_posns = contig_map.second;
            if (contig_set.find(contig_id) == contig_set.end()) {
                continue;
            }

            //initialize alignment data structures
            vector<float> contig_ncp_vector = non_collapse_prob_map_contig[contig_id];
            vector<vector<float>> S(contig_posns.size(),vector<float>(x_posns.size()-1));
            vector<vector<array<int,2>>> previous(contig_posns.size(),vector<array<int,2>>(x_posns.size()-1));

            if (!local_aln) { //handles standard and tip
                init_semiglobal_aln(S,x_posns.size(),contig_posns.size(),x);

            } else {
                init_semiglobal_aln(S,x_posns.size(),contig_posns.size(),x);
            }

            set<int> temp;
            dp_align(S, previous, contig_posns, x_posns, false, lookback, temp, seg_ncp_vector, contig_ncp_vector, swap_b_x);
            array<int, 2> best_start;
            if (!local_aln)
                best_start = sg_aln_backtrack_start(S,contig_posns.size(),x_posns.size());
            else if (tip_aln) {
                best_start = sg_aln_backtrack_start(S,contig_posns.size(),x_posns.size());
            } else {
                best_start = local_backtrack_start(S, contig_posns.size(), x_posns.size());
            }

            float curr_best_score = S[best_start[0]][best_start[1]];
            if (detection) {
                //determines how many alns to get for this reference sequence
                int alns_to_get = int(roundf(n_detect_scores * float(contig_posns.size()) / float(total_labels))) + 1;
                multi_backtrack_scores(S, previous, x, contig_id, contig_posns.size(), x_posns.size(), alns_to_get,
                                       scoring_vector);
            } else {
                scoring_vector.emplace_back(x, contig_id, curr_best_score);
            }
        }

    }
    return scoring_vector;
}

void run_SA_fitting(map<int,vector<float>> cmaps_segs, map<int,vector<float>> cmaps_contigs,
        map<int,vector<float>> non_collapse_prob_map, map<int,vector<float>> non_collapse_prob_map_contig,
        vector<pair<int,int>> pairs) {

//    cmap_map_to_string(cmaps_segs);
//    cout << "\n";
//    cmap_map_to_string(cmaps_contigs);

    map<pair<int,int>,set<int>> curr_aligned_labs;
    map<int,set<int>> discovered_contig_used_labels;
    for (auto p: pairs) {
        int x = p.first;
        int contig_id = p.second;
        curr_aligned_labs[make_pair(x,contig_id)] = set<int>();

        vector<float> x_posns = cmaps_segs[x];
        vector<float> contig_posns = cmaps_contigs[contig_id];

        vector<float> seg_ncp_vector = non_collapse_prob_map[x];
        vector<float> contig_ncp_vector = non_collapse_prob_map_contig[contig_id];

        vector<vector<float>> S(contig_posns.size(), vector<float>(x_posns.size() - 1));
        vector<vector<array<int, 2>>> previous(contig_posns.size(), vector<array<int, 2>>(x_posns.size() - 1));

        init_fitting_aln(S, x_posns.size(),contig_posns.size(),x);
        dp_align(S, previous, contig_posns, x_posns, true, lookback, curr_aligned_labs[(make_pair(x,contig_id))], seg_ncp_vector,
                contig_ncp_vector, swap_b_x);

        int contig_end = contig_posns.size()- 2;
        int seg_end = x_posns.size()-2;
        array<int, 2> best_start = {contig_end,seg_end};
        vector<tuple<int, int, float>> aln_list = get_aln_list(S, previous, best_start);

        float curr_best_score = S[best_start[0]][best_start[1]];
        string seg_id = to_string(abs(x));
        string outname = sample_prefix + "_" + to_string(contig_id) + "_" + seg_id + "_fitting_aln.txt";
        print_alignment(S, previous, aln_list, contig_id, x, cmaps_segs, outname, discovered_contig_used_labels, curr_best_score);

    }
}

/*Method to do the alignments
* Takes a list of the contigs and segments to align, writes the stuff to files on its own
*/
map<int,set<int>> run_SA_aln(map<int,vector<float>> cmaps_segs, map<int,vector<float>> cmaps_contigs,
                map<int,vector<float>> non_collapse_prob_map, map<int,vector<float>> non_collapse_prob_map_contig,
                set<pair<int,int>> seg_contig_pairs, map<int,set<int>> contig_used_label_map, map<int,float> score_thresholds) {


    map<int,set<int>> discovered_contig_used_labels;
    map<pair<int,int>,set<int>> curr_aligned_labs;

    for (auto e: cmaps_contigs) {
        discovered_contig_used_labels[e.first] = set<int>();
    }

    //create variables for the segment and contig iterations
    for (auto p: seg_contig_pairs) {
        int x = p.first;
        int contig_id = p.second;
        curr_aligned_labs[make_pair(x,contig_id)] = set<int>();
        if (tip_aln) {
            curr_aligned_labs[make_pair(x,contig_id)] = contig_used_label_map[contig_id];
        }

        int min_aln_labels = 6;
        vector<float> x_posns = cmaps_segs[x];
        vector<float> contig_posns = cmaps_contigs[contig_id];
        vector<float> seg_ncp_vector = non_collapse_prob_map[x];
        vector<float> contig_ncp_vector = non_collapse_prob_map_contig[contig_id];

        //initialize alignment data structures
        float curr_best_score = score_thresholds[x];
        float exp_thresh = score_thresholds[x];
        int x_aligned_count = 0;
        float mean_thresh,med_thresh;
        if (tip_aln) {
            med_thresh = 8500;
            mean_thresh = 8750;
        } else {
            med_thresh = 8000;
            mean_thresh = 7000;
        }

//        if ((abs(x) == 2 && contig_id == 80) || (x == -48 && contig_id == 80)) {
//            cout << x << " " << contig_id << " " << curr_best_score << " " << "\n";
//        }

        //Run the alignment while the score exceeds the threshold
        while (curr_best_score >= exp_thresh && x_aligned_count < max_seg_contig_alns) {
            vector<vector<float>> S(contig_posns.size(), vector<float>(x_posns.size() - 1));
            vector<vector<array<int, 2>>> previous(contig_posns.size(), vector<array<int, 2>>(x_posns.size() - 1));
            if (!local_aln) {
                init_semiglobal_aln(S,x_posns.size(),contig_posns.size(),x);

            } else {
                init_local_aln(S,x_posns.size(),contig_posns.size(),x);
            }

            dp_align(S, previous, contig_posns, x_posns, false, lookback, curr_aligned_labs[(make_pair(x,contig_id))], seg_ncp_vector,
                     contig_ncp_vector, swap_b_x);

            array<int, 2> best_start;
            if (!local_aln) {
                best_start = sg_aln_backtrack_start(S, contig_posns.size(), x_posns.size());
            } else {
                best_start = local_backtrack_start(S, contig_posns.size(), x_posns.size());
            }
            vector<tuple<int, int, float>> aln_list = get_aln_list(S, previous, best_start);
            float exp_score_thresh_lab, exp_score_thresh_len;
            tie(exp_score_thresh_lab,exp_score_thresh_len) = compute_partial_score_threshold(score_thresholds[x], aln_list, cmaps_contigs[contig_id]);
            exp_thresh = fmax(exp_score_thresh_lab,exp_score_thresh_len);
            curr_best_score = S[best_start[0]][best_start[1]];

            string tip_aln_status;
            //if it passes the e-value write it, and do not update the contig
            if (curr_best_score > exp_thresh) {
                x_aligned_count += 1;

                for (auto &e: aln_list) {
                    curr_aligned_labs[(make_pair(x,contig_id))].insert(get<0>(e));
                }

                //make sure it's actually a tip, not a middle.
                if (tip_aln) {
                    tip_aln_status = "_tip";
                    if (((contig_posns.size()-2) - get<0>(aln_list[0]) > 3) && get<0>(aln_list.back()) > 2) {
                        continue;
                    }
                }

                if (aln_list.size() < 3 || (!tip_aln && aln_list.size() < min_aln_labels)) {
                    continue;
                }

                //mean and median checks
                float mean = get<2>(aln_list[0])/(aln_list.size()-1);
                float median = compute_aln_median(aln_list);

                if (mean < mean_thresh || median < med_thresh) {
                    continue;

                } else if (aln_list.size() == 3) {
                    float diff1 = (get<2>(aln_list[0]) - get<2>(aln_list[1]));
                    float diff2 = (get<2>(aln_list[1]) - get<2>(aln_list[2]));
                    if (diff1 < med_thresh || diff2 < med_thresh) {
                        continue;
                    }
                }

                string seg_id = to_string(abs(x));
                if (x < 0) {
                    seg_id += "_r";
                }
                string outname = sample_prefix + "_" + to_string(contig_id) + "_" + seg_id + "_" +
                                 to_string(x_aligned_count) + flipped + tip_aln_status + "_aln.txt";
                print_alignment(S, previous, aln_list, contig_id, x, cmaps_segs, outname, discovered_contig_used_labels, curr_best_score);
            }
        }
    }
    return discovered_contig_used_labels;
}

 /*
  * method to parse the args
  */
 tuple<string,bool,bool> parse_args(int argc, char *argv[]) {
     //Parse command line arguments
     string contig_list_file;
     string used_labels_file;
     bool limit_lookback = true;
     bool do_tip = true;
     //-------------------------------------------------------
     for (int i = 3; i < argc; ++i) {
         if (string(argv[i]) == "-nl") {
             limit_lookback = false; //from SegAligner.h

         } else if (string(argv[i]) == "-no_tip_aln") {
            do_tip = false;

            //TODO: REDUNDANT ARGS BELOW
         } else if (string(argv[i]) == "-detection") {
             local_aln = true;
             detection = true;
             swap_b_x = true; //from SegAligner.h
             do_tip = false;
             detection_label = "_detection";
             p_val = p_val_RG;
             cout << "detection mode ON, local ON, swap_b_x ON" << endl;

//         } else if (string(argv[i]) == "-alnref") {
//             local_aln = true;
//             swap_b_x = true; //from SegAligner.h
//             flipped = "_flipped";
//             do_tip = false;
//             cout << "detection mode ON (local ON)" << endl;

         } else if (string(argv[i]).rfind("-prefix=", 0) == 0) {
             sample_prefix = string(argv[i]).substr(string(argv[i]).find('=') + 1);

         } else if (string(argv[i]).rfind("-contig_list=", 0) == 0) {
             contig_list_file = string(argv[i]).substr(string(argv[i]).find('=') + 1);

         } else if (string(argv[i]) == "-local") {
             local_aln = true;
             cout << "Local alignment mode" << endl;

         } else if (string(argv[i]) == "-fitting") {
             fitting_aln = true;
             cout << "Fitting alignment mode. Multithreading will be off." << endl;
             cout << "" << endl;

         } else if (string(argv[i]).rfind("-nthreads=", 0) == 0) {
             n_threads = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

         } else if (string(argv[i]).rfind("-min_labs=", 0) == 0) {
             min_map_len = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

         } else if (string(argv[i]).rfind("-n_detect_scores=", 0) == 0) {
             //need num scores if detection is true
             n_detect_scores = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));

         } else if (string(argv[i]).rfind("-gen=", 0) == 0) {
             inst_gen = stoi(string(argv[i]).substr(string(argv[i]).find('=') + 1));
             if (inst_gen != 1 && inst_gen != 2) {
                 cout << "1 or 2 are the only valid arguments for -gen";
                 exit(1);
             }
         }

         if (!flipped.empty()) {
             do_tip = false;
         }
     }
     //TODO: Print a summary of the positional arguments given
    return make_tuple(contig_list_file,limit_lookback,do_tip);

 }

 set<int> get_contig_set(string contig_list_file,map<int,vector<float>> cmaps_contigs) {
     set<int> contig_set;
     if (!contig_list_file.empty()) {
         contig_set = parse_contig_list(contig_list_file);
     } else {
         for (auto &key: cmaps_contigs) {
             contig_set.insert(key.first);
         }
     }
     return contig_set;
 }


/*
 *MAIN
 */
int main (int argc, char *argv[]) {
    if (argc < 2) {
        cout << "wrong number of arguments" << endl;
        exit(1);
    }

    //timing
    clock_t begin = clock();
    double elapsed_secs;
    clock_t end;

    //positional arguments
    string ref_cmap_file = argv[1];
    string contig_cmap_file = argv[2];

    //Parse command line args
    string contig_list_file;
    bool limit_lookback, do_tip;
    tie(contig_list_file, limit_lookback, do_tip) = parse_args(argc,argv);

    if (!limit_lookback) {
        lookback = 9999999;
        cout << "Alignment banding OFF. \n";
    }

    //make segs cmap
    map<int,vector<float>> cmaps_segs_raw;

    parse_cmap(ref_cmap_file,cmaps_segs_raw);
    //cmap_map_to_string(cmap_map) //DEBUG
    //add reverse segs
    map<int, vector<float>>cmaps_segs;
    if (!fitting_aln) {
        cmaps_segs = make_reverse_cmap(cmaps_segs_raw, min_map_len);
    } else {
        cmaps_segs = cmaps_segs_raw;
    }

    //make contigs cmaps
    map<int,vector<float>> cmaps_contigs;
    parse_cmap(contig_cmap_file, cmaps_contigs);

    //find collapse probs
    map<int,vector<float>> non_collapse_prob_map = non_collapse_probs(cmaps_segs,inst_gen);
    map<int,vector<float>> non_collapse_prob_map_contig = non_collapse_probs(cmaps_contigs,inst_gen);

    //-------------------------------------------------------

    //[Alignment] Get contigs to align with
    set<int> contig_set = get_contig_set(contig_list_file, cmaps_contigs);

    cout << "Running SegAligner on " << cmaps_segs.size() << " segments and " << cmaps_contigs.size() << " contigs.\n";

    //if fitting alignment mode, just align, report and quit.
    if (fitting_aln) {
        n_threads = 1;
        vector<pair<int,int>> pairs;
        for (auto e: contig_set) {
            for (auto x: cmaps_segs) {
                pairs.emplace_back(x.first,e);
            }
        }

        run_SA_fitting(cmaps_segs, cmaps_contigs, non_collapse_prob_map, non_collapse_prob_map_contig,pairs);
        end = clock();
        elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        printf("elapsed seconds for run %f\n",elapsed_secs);
        return 0;
    }

    //Get scoring distribution
    cout << "Computing scoring thresholds.\n";
    vector<future<vector<tuple<int,int,float>>>> future_scores;
    vector<set<int>> chunked_contig_sets(n_threads,set<int>());
    //chunk the relevant contigs per thread
    int count = 0;
    for (auto e: contig_set) {
        chunked_contig_sets[(count % n_threads)].insert(e);
        count++;
    }

    //[RUN SCORING]
    for (int i = 0; i < n_threads; i++) {
        future_scores.push_back(async(launch::async, run_SA_score, cmaps_segs, cmaps_contigs, non_collapse_prob_map,
                non_collapse_prob_map_contig, chunked_contig_sets[i]));
    }


    //Gather results
    vector<tuple<int,int,float>>full_scoring_vector;
    for (auto &e: future_scores) {
        vector<tuple<int,int,float>> temp = e.get();
        full_scoring_vector.insert(full_scoring_vector.end(),temp.begin(),temp.end());
    }

    write_all_scores(full_scoring_vector, sample_prefix);

    map<int,float>score_thresholds = compute_score_thresholds(full_scoring_vector,cmaps_segs,cmaps_contigs,p_val);
    map<int,float>tip_thresholds = compute_score_thresholds(full_scoring_vector,cmaps_segs,cmaps_contigs,p_val_tip);
    write_score_thresholds(score_thresholds, sample_prefix + detection_label);
    cout << "Completed generating scores.\n";
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("elapsed seconds for scoring %f\n",elapsed_secs);

    //[PREP ALIGNMENTS]
    vector<pair<int,int>> pairs_list;
    //which pairs exceeded scoring dist
    for (auto e: full_scoring_vector) {
        int x = get<0>(e);
        float curr_score = get<2>(e);
        if (curr_score >= score_thresholds[x]) {
            pairs_list.emplace_back(get<0>(e),get<1>(e));
        }
    }

    //chunk the pairings passing threshold
    vector<set<pair<int,int>>> chunked_align_sets(n_threads,set<pair<int,int>>());
    //chunk the relevant contigs per thread
    for (int i = 0; i < pairs_list.size(); i++) {
        chunked_align_sets[(i % n_threads)].insert(pairs_list[i]);
    }

    //[For Tip alignment] Get the set of labels for each contig which have already been aligned
    map<int,set<int>> contig_used_label_map;
    for (auto &e: cmaps_contigs) {
        contig_used_label_map[e.first] = set<int>();
    }

    //[RUN ALIGNMENTS]
    cout << "Performing alignments\n";
    vector<future<map<int,set<int>>>> futs;
    map<int,set<int>> temp;
    for (int i = 0; i < n_threads; i++) {
        futs.push_back(async(launch::async, run_SA_aln, cmaps_segs, cmaps_contigs, non_collapse_prob_map,
                non_collapse_prob_map_contig, chunked_align_sets[i], temp, score_thresholds));
    }

    //TODO: Switch to get and merge map
    for (auto &f: futs) {
        temp = f.get();
        for (auto &x: temp) {
            contig_used_label_map[x.first].insert(x.second.begin(),x.second.end());
        }
    }

    cout << "Finished standard alignment. \n";

    //[TIP ALIGNMENT]
    //Score then align tips
    if (do_tip) {
        cout << "Doing tip alignment.\n";

        set<pair<int,int>> tip_pairs_set;
        tip_aln = true;
//        local_aln = true;
        //Get contigs with alns and pair with all segs
        for (auto e: pairs_list) {
            //don't make the pair if there weren't any alignments to the contigs. check the used label set
            if (!contig_used_label_map[e.second].empty()) {
//                for (auto x: pairs_list) {
                for (auto x: cmaps_segs) {
                //adds the seg_id and contig_id so all relevant segs are paired with all relevant contigs
                    tip_pairs_set.insert(make_pair(x.first, e.second));
                    tip_pairs_set.insert(make_pair(-x.first, e.second));
                }
            }
        }

        //chunk the tip pairs
        vector<set<pair<int,int>>> chunked_tip_pairs(n_threads,set<pair<int,int>>());
        count = 0;
        for (auto e: tip_pairs_set) {
            chunked_tip_pairs[(count % n_threads)].insert(e);
            count++;
        }


        //Run SA align
        cout << "Performing alignments\n";
        vector<future<map<int,set<int>>>> tip_futs;
        for (int i = 0; i < n_threads; i++) {
            tip_futs.push_back(async(launch::async, run_SA_aln, cmaps_segs, cmaps_contigs, non_collapse_prob_map,
                                 non_collapse_prob_map_contig, chunked_tip_pairs[i], contig_used_label_map, tip_thresholds));

        }
        for (auto &f: tip_futs) {
            f.wait();
        }
        cout << "Finished tip alignments.\n" << endl;

    }

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("elapsed seconds for run %f\n",elapsed_secs);

    return 0;
}