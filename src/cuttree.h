#pragma once

#include <unordered_map>

namespace cuttree {

// clusters are of edges, so size = 2 => 3 nodes
constexpr int MIN_CLUSTER_SIZE = 2;

struct EdgeComponent
{
    double d;
    int from, to, cluster_num;
};

struct TreeDat
{
    std::vector <EdgeComponent> edges;
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

void fill_edges (TreeDat &tree,
        const std::vector <int> &from,
        const std::vector <int> &to,
        Rcpp::NumericVector &d);
double calc_ss (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
double calc_covsum (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
size_t cluster_size (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
std::unordered_set <int> build_one_tree (std::vector <EdgeComponent> &edges);

TwoSS sum_component_ss (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree_edges, const bool shortest);
BestCut find_min_cut (const TreeDat &tree, const int cluster_num,
        const bool shortest);

} // end namespace cuttree

Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl,
        const bool shortest);
