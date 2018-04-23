#pragma once

// --------- AVERAGE LINKAGE CLUSTER ----------------

#include "bst.h"

#include <unordered_set>

/* The main matrices (contig, num_edges, dmat, avg_dist) are all referenced by
 * direct indices throughout, not by vertex numbers. The latter are mapped to
 * the former by vert2index_map. Note that index2vert_map is not used for this
 * routine, but exists as dummy to pass to `sets_init`
 *
 * The index2cl and cl2index then associate those indices with clusters which
 * are themselves also direct indices into the matrices. Cluster merging simply
 * re-directs multiple indices onto the same cluster (index) numbers.
 *
 * The binary tree only returns minimal distances which need to be associated
 * with particular pairs of clusters. This is done with the final map,
 * edgewt2idx_pair, where the pair of indices is into clusters, requiring this
 * map to be constantly updated. This updating requires in turn a reverse map,
 * idx2edgewt, so that the weight associated with any pre-merge cluster can
 * be obtained, and the edgewt2idx clusters for that weight updated.
 */
struct ALKDat
{
    unsigned int n;

    std::unordered_map <double,
        std::pair <unsigned int, unsigned int> > edgewt2idx_pair_map;
    std::unordered_map <unsigned int, std::unordered_set <double> >
        idx2edgewt_map; // all wts associated with that cluster

    arma::Mat <unsigned short> contig_mat, num_edges;
    arma::Mat <double> dmat, avg_dist;

    uint_map_t index2cl_map, vert2index_map, index2vert_map;
    uint_set_map_t cl2index_map;
};

void alk_init (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void update_edgewt_maps (ALKDat &alk_dat,
        unsigned int l, unsigned int m);

int alk_step (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
