#pragma once

#include <unordered_map>

constexpr int MIN_CLUSTER_SIZE = 3;

struct EdgeComponent
{
    double d;
    int from, to, cluster_num;
};

struct BestCut
{
    int pos, n1, n2;
    double ss_diff, ss1, ss2;
    std::unordered_set <int> nodes;
};

struct TwoSS // 2 sums-of-squares values
{
    double ss1, ss2;
    int n1, n2; // sizes of clusters
};

void fill_edges (std::vector <EdgeComponent> &edges,
    const std::vector <std::string> &from,
    const std::vector <std::string> &to,
    const Rcpp::NumericVector &d);
double calc_ss (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
size_t cluster_size (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
std::unordered_set <int> build_one_tree (std::vector <EdgeComponent> &edges);
TwoSS sum_component_ss (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree);
BestCut find_min_cut (const std::vector <EdgeComponent> &edges,
        const int cluster_num);

Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl);
