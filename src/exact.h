#pragma once

// --------- EXACT CLUSTER ----------------

struct EXDat
{
    unsigned int n;

    std::vector <oneEdge> edges_all, edges_nn;

    arma::Mat <unsigned short> contig_mat;
    arma::Mat <double> dmax;

    uint_map_t index2cl_map, vert2index_map, index2vert_map;
    uint_set_map_t cl2index_map;
};

void clexact_init (EXDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

unsigned int clexact_step (EXDat &clexact_dat, unsigned int i);

Rcpp::IntegerVector rcpp_exact (
        const Rcpp::DataFrame gr);
