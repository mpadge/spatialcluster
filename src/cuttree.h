#pragma once

#include <unordered_map>

struct EdgeComponent
{
    double d;
    int from, to, cluster_num;
};

struct BestCut
{
    size_t n1, n2; // TODO: Delete those
    int pos;
    double variance;
};

double simple_variance (const std::vector <double> &x);
double calc_variance (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
size_t cluster_size (const std::vector <EdgeComponent> &edges,
        const int cluster_num);
std::unordered_set <int> build_one_tree (std::vector <EdgeComponent> &edges);
double sum_component_variances (const std::vector <EdgeComponent> &edges,
        const std::unordered_set <int> &tree);
BestCut find_min_cut (std::vector <EdgeComponent> &edges,
        const int cluster_num);

Rcpp::IntegerVector rcpp_cut_tree (const Rcpp::DataFrame tree, const int ncl);
