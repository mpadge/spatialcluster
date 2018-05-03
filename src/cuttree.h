#pragma once

#include <unordered_map>

struct EdgeComponent
{
    double d;
    int from, to, cluster_num;
};

struct BestCut
{
    int pos;
    double ss1, ss2;
};

struct TwoSS // 2 sums-of-squares values
{
    double ss1, ss2;
};

double calc_ss (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
size_t cluster_size (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
std::unordered_set <int> build_one_tree (std::vector <EdgeComponent> &edges);
TwoSS sum_component_ss (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree);
BestCut find_min_cut (std::vector <EdgeComponent> &edges,
        const int cluster_num);

Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl);
