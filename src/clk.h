#pragma once

#include "utils.h"

// --------- COMPLETE LINKAGE CLUSTER ----------------

namespace clk {

struct CLKDat {
    bool shortest;
    size_t n;

    std::vector <utils::OneEdge> edges_all, edges_nn;

    arma::Mat <int> contig_mat;
    arma::Mat <double> dmat;

    int2indxset_map_t cl2index_map;

    indx2int_map_t index2cl_map, index2vert_map;
    int2indx_map_t vert2index_map;
};

void clk_init (CLKDat &clk_dat,
        Rcpp::IntegerVector from_full,
        Rcpp::IntegerVector to_full,
        Rcpp::NumericVector d_full,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

size_t clk_step (CLKDat &clk_dat, size_t i);

} // end namespace clk

Rcpp::IntegerVector rcpp_clk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr,
        const bool shortest,
        const bool quiet);
