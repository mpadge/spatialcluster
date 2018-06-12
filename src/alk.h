#pragma once

// --------- AVERAGE LINKAGE CLUSTER ----------------

#include "bst.h"

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

namespace alk {

struct ALKDat
{
    size_t n;

    std::unordered_map <double,
        std::pair <index_t, index_t> > edgewt2idx_pair_map;
    std::unordered_map <index_t, std::unordered_set <double> >
        idx2edgewt_map; // all wts associated with that cluster

    arma::Mat <int> contig_mat, num_edges;
    arma::Mat <double> dmat, avg_dist;

    int2indxset_map_t cl2index_map;
    indx2int_map_t index2cl_map, index2vert_map;
    int2indx_map_t vert2index_map;
};

void alk_init (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void update_edgewt_maps (ALKDat &alk_dat, index_t l, index_t m);

size_t alk_step (ALKDat &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d,
        bool distances);

} // end namespace alk

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
