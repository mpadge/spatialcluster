#pragma once

int find_shortest_connection (
        Rcpp::IntegerVector &from,
        Rcpp::IntegerVector &to,
        Rcpp::NumericVector &d,
        arma::Mat <double> &d_mat,
        uint_set_map_t &cl2edge_map,
        int cfrom,
        int cto);

void merge_clusters (
        arma::Mat <unsigned short> &contig_mat,
        uint_map_t &vert2cl_map,
        uint_set_map_t &cl2vert_map,
        int merge_from,
        int merge_to);

Rcpp::IntegerVector rcpp_slk (
        const Rcpp::DataFrame gr_full,
        const Rcpp::DataFrame gr);
