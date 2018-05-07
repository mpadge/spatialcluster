#pragma once

#include "utils.h"

// --------- EXACT MERGE ----------------
//
// Merge the clusters generated by the rcpp_exact_initial. Seperate class and
// routines to allow results from rcpp_exact_initial to be returned and cached
// for subsequent re-merging.

struct OneCluster
{
    int id, n;
    double dist_sum, dist_max;
    std::vector <OneEdge> edges;
};

struct OneMerge
{
    int cli, clj;
    double merge_dist;
};

struct EXMerge
{
    // cl2index_map is from cluster numbers to indices in clusters
    std::unordered_map <int, unsigned int> cl2index_map;
    std::vector <OneCluster> clusters;
    std::vector <OneEdge> edges; // edges between clusters
    std::vector <OneMerge> merges;
};

void rcpp_exmerge_init (EXMerge &cldat);
OneMerge rcpp_exmerge_merge (EXMerge &cldat,
        int clfrom_i,
        int clto_i,
        unsigned int ei);
void rcpp_exmerge_single (EXMerge &cldat);
void rcpp_exmerge_avg (EXMerge &cldat);
void rcpp_exmerge_max (EXMerge &cldat);

Rcpp::IntegerVector rcpp_exact_merge (
        const Rcpp::DataFrame gr,
        const int ncl,
        const std::string method);
