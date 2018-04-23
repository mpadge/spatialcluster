#pragma once

// --------- AVERAGE LINKAGE CLUSTER ----------------

#include "bst.h"

#include <unordered_set>

/* This relies on the CLDAT structure in common.h, with index2vert_map not used
 * for this routine, although it exists as dummy to pass to `sets_init`
 *
 * The binary tree only returns minimal distances which need to be associated
 * with particular pairs of clusters. This is done with the final map,
 * edgewt2idx_pair, where the pair of indices is into clusters, requiring this
 * map to be constantly updated. This updating requires in turn a reverse map,
 * idx2edgewt, so that the weight associated with any pre-merge cluster can
 * be obtained, and the edgewt2idx clusters for that weight updated.
 */

void alk_init (CLDAT &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

void update_edgewt_maps (CLDAT &alk_dat,
        unsigned int l, unsigned int m);

int alk_step (CLDAT &alk_dat,
        BinarySearchTree &tree,
        Rcpp::IntegerVector from,
        Rcpp::IntegerVector to,
        Rcpp::NumericVector d);

Rcpp::IntegerVector rcpp_alk (
        const Rcpp::DataFrame gr);
