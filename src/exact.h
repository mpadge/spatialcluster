#pragma once

// --------- EXACT CLUSTER ----------------

/* All `index2cl` values are initally set to -1, and there are no `cl2index`
 * values.  There is also a single binary vector of `index_in_cluster`, initialy
 * set to `false`.
 */
struct EXDat
{
    unsigned int n;

    std::vector <oneEdge> edges; // nearest neighbour edges only
    std::vector <bool> index_in_cluster;

    uint_map_t index2cl_map, vert2index_map, index2vert_map;
    uint_set_map_t cl2index_map;
};

void clexact_init (EXDat &clexact_dat,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void assign_first_edge (EXDat &clexact_dat);

unsigned int clexact_step (EXDat &clexact_dat, unsigned int ei,
        unsigned int clnum);

Rcpp::IntegerVector rcpp_exact (
        const Rcpp::DataFrame gr);
