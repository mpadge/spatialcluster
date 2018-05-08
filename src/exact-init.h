#pragma once

#include "utils.h"

// --------- EXACT CLUSTER ----------------

/* All `index2cl` values are initally set to -1, and there are no `cl2index`
 * values.  There is also a single binary vector of `index_in_cluster`, initialy
 * set to `false`.
 */

namespace ex_init {

struct ExInitDat
{
    size_t n;

    std::vector <utils::OneEdge> edges; // nearest neighbour edges only
    std::vector <bool> index_in_cluster;

    int2int_map_t vert2cl_map;
    int2indx_map_t vert2index_map;
    indx2int_map_t index2cl_map, index2vert_map;
    int2intset_map_t cl2index_map;
};

void init (ExInitDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void assign_first_edge (ExInitDat &clexact_dat);

int step (ExInitDat &clexact_dat, const index_t ei,
        const int clnum);

void fill_cl_edges (ExInitDat &clexact_dat, arma::Mat <double> &cl_edges,
        int num_clusters);

} // end namespace ex_init

Rcpp::IntegerVector rcpp_exact_initial (
        const Rcpp::DataFrame gr);
